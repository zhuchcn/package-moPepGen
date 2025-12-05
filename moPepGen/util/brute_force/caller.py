""" Refactored BruteForceVariantPeptideCaller with composition """
import copy
import math
from typing import Iterable, List, Set, Tuple
from itertools import combinations
from Bio import SeqUtils
from Bio.Seq import Seq
from moPepGen import gtf, seqvar, aa, dna, params, constant
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool
from moPepGen.seqvar.VariantRecordWithCoordinate import VariantRecordWithCoordinate
from .orf_tracker import ORFTracker
from .sequence_builder import VariantSequenceBuilder
from .effect_analyzer import VariantEffectAnalyzer
from .validator import VariantCompatibilityValidator


class BruteForceVariantPeptideCaller:
    """
    Variant peptide caller using the brute force algorithm.
    
    This class orchestrates peptide calling by delegating to helper classes:
    - VariantSequenceBuilder: Generates variant sequences
    - VariantEffectAnalyzer: Analyzes variant effects
    - VariantCompatibilityValidator: Validates variant combinations
    """
    
    def __init__(self, reference_data:params.ReferenceData=None,
            cleavage_params:params.CleavageParams=None,
            variant_pool:VariantRecordPool=None,
            canonical_peptides=None, tx_id:str=None,
            tx_model:gtf.TranscriptAnnotationModel=None,
            tx_seq:dna.DNASeqRecordWithCoordinates=None, gene_seq:dna.DNASeqRecord=None,
            start_index:int=None, variant_peptides:Set[str]=None,
            w2f:bool=False, selenocysteine_termination:bool=False):
        """ Constructor """
        self.reference_data = reference_data
        self.cleavage_params = cleavage_params
        self.variant_pool = variant_pool or VariantRecordPool()
        self.canonical_peptides = canonical_peptides
        self.tx_id = tx_id
        self.tx_model = tx_model
        self.tx_seq = tx_seq
        self.gene_seq = gene_seq
        self.start_index = start_index
        self.variant_peptides = variant_peptides or set()
        self.selenocysteine_termination = selenocysteine_termination
        self.w2f = w2f
        
        # Helper classes (lazy initialization)
        self._sequence_builder = None
        self._effect_analyzer = None
        self._validator = None
    
    @property
    def sequence_builder(self) -> VariantSequenceBuilder:
        """ Lazy initialization of sequence builder """
        if self._sequence_builder is None:
            self._sequence_builder = VariantSequenceBuilder(
                self.reference_data, self.tx_model, self.tx_seq, 
                self.variant_pool, self.tx_id
            )
        return self._sequence_builder
    
    @property
    def effect_analyzer(self) -> VariantEffectAnalyzer:
        """ Lazy initialization of effect analyzer """
        if self._effect_analyzer is None:
            self._effect_analyzer = VariantEffectAnalyzer(
                self.tx_seq, self.tx_model, self.tx_id
            )
        return self._effect_analyzer
    
    @property
    def validator(self) -> VariantCompatibilityValidator:
        """ Lazy initialization of validator """
        if self._validator is None:
            self._validator = VariantCompatibilityValidator(
                self.reference_data, self.tx_model, self.tx_seq,
                self.variant_pool, self.tx_id
            )
        return self._validator
    
    def create_canonical_peptide_pool(self):
        """ Create canonical peptide pool. """
        proteome = self.reference_data.proteome
        par = self.cleavage_params
        self.canonical_peptides = proteome.create_unique_peptide_pool(
            anno=self.reference_data.anno,
            rule=par.enzyme,
            exception=par.exception,
            miscleavage=par.miscleavage,
            min_mw=par.min_mw,
            min_length=par.min_length,
            max_length=par.max_length
        )
    
    def get_start_index(self):
        """ Get the "start index" used for filtering variants. """
        if self.tx_seq.orf:
            self.start_index = self.tx_seq.orf.start + 3
        else:
            self.start_index = 3
    
    def peptide_is_valid(self, peptide:str, denylist:List[str], check_canonical) -> bool:
        """ Check whether the peptide is valid """
        if check_canonical \
                and self.canonical_peptides \
                and peptide in self.canonical_peptides:
            return False
        if peptide in denylist:
            return False
        min_len = self.cleavage_params.min_length
        max_len = self.cleavage_params.max_length
        min_mw = self.cleavage_params.min_mw
        return min_len <= len(peptide) <= max_len \
            and SeqUtils.molecular_weight(peptide, 'protein') >= min_mw
    
    def load_relevant_variants(self, pool:VariantRecordPool):
        """ Load relevant variants. """
        if self.tx_id not in pool:
            return
        
        for variant in pool[self.tx_id].transcriptional:
            self.variant_pool.add_transcriptional_variant(variant)
        for variant in pool[self.tx_id].intronic:
            self.variant_pool.add_intronic_variant(variant)
        for variant in pool[self.tx_id].fusion:
            self.variant_pool.add_fusion_variant(variant)
            right_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
            if right_tx_id not in pool:
                continue
            for v in pool[right_tx_id].transcriptional:
                self.variant_pool.add_transcriptional_variant(v)
            for v in pool[right_tx_id].intronic:
                self.variant_pool.add_intronic_variant(v)
        for variant in pool[self.tx_id].circ_rna:
            self.variant_pool.add_circ_rna(variant)
    
    def get_variants(self, tx_id:str, start:int, end:int,
            variant_ids:Iterable[str]=None) -> List[VariantRecord]:
        """ Load variant records associated with the particular transcript. """
        series = self.variant_pool[tx_id]
        variants = []
        for variant in series.transcriptional:
            if variant.location.start < start -1:
                continue
            if self.tx_model.is_mrna_end_nf() and variant.location.end <= end:
                continue
            if variant_ids:
                if variant.id not in variant_ids:
                    continue
            if variant.location.start == start - 1:
                variant.to_end_inclusion(self.tx_seq)
            variants.append(variant)
        return variants
    
    @staticmethod
    def find_prev_cds_start_same_frame(cds_start:int, cds_start_positions:List[int]):
        """ find the previous cds start site in the same reading frame. """
        if cds_start == 0:
            return -1
        reading_frame_index = cds_start % 3
        for site in cds_start_positions[::-1]:
            if site >= cds_start:
                continue
            if site % 3 == reading_frame_index:
                return site
        return -1
    def get_effective_variants(self, lhs:int, rhs:int,
            lrange:Tuple[int,int], rrange:Tuple[int,int], cds_start:int,
            variants:List[seqvar.VariantRecordWithCoordinate],
            variants_stop_lost:List[Tuple[bool,bool,bool]],
            variants_stop_gain:List[Tuple[bool,bool,bool]],
            variants_silent_mutation:List[Tuple[bool,bool,bool]]
            ) -> List[VariantRecordWithCoordinate]:
        """ Check whether the given range of the transcript has any variant
        associated. """
        effective_variants:List[seqvar.VariantRecordWithCoordinate] = []
        offset = 0
        query = FeatureLocation(start=lhs, end=rhs)
        start_loc = FeatureLocation(start=cds_start, end=cds_start + 3)
        rf_index = cds_start % 3
        frames_shifted = 0
        upstream_indels:List[VariantRecordWithCoordinate] = []
        is_coding = self.tx_model.is_protein_coding \
            and not any(v.variant.is_circ_rna() for v in variants)
        for i, variant_coordinate in enumerate(variants):
            variant = variant_coordinate.variant
            loc = variant_coordinate.location
            if variant.type == 'Insertion':
                loc = FeatureLocation(start=loc.start+1, end=loc.end)
            if loc.start > rhs + 3:
                break
            is_start_gain = start_loc.overlaps(loc)
            is_frameshifting = cds_start < loc.start < lhs and variant.is_frameshifting()
            if cds_start < loc.start < lhs:
                frames_shifted = (frames_shifted + variant.frames_shifted()) % 3
                if variant.is_frameshifting() \
                        and not (variant.is_fusion() or variant.is_circ_rna()):
                    upstream_indels.append(variant_coordinate)
            is_cleavage_gain = (
                loc.overlaps(FeatureLocation(start=lrange[0], end=lrange[1]))
                or loc.overlaps(FeatureLocation(start=rrange[0], end=rrange[1]))
            ) \
                if cds_start != lhs \
                else loc.overlaps(FeatureLocation(start=rhs, end=rhs + 3))

            is_stop_lost = variants_stop_lost[i][cds_start % 3] \
                and cds_start < variants[i].location.start < rhs

            is_silent_mutation = variants_silent_mutation[i][cds_start % 3] \
                and variants[i].location.start > cds_start

            is_stop_gain = variants_stop_gain[i][cds_start % 3] \
                and int((variants[i].location.start - rf_index) / 3) \
                    == int((rhs - rf_index)/3)

            if (loc.overlaps(query)
                        or is_start_gain
                        or (is_frameshifting and not is_coding)
                        or is_cleavage_gain
                        or is_stop_lost
                        or is_stop_gain ) \
                    and not is_silent_mutation:
                effective_variants.append(variant_coordinate)
            offset += len(variant.alt) - len(variant.ref)

        if frames_shifted > 0:
            for v in upstream_indels:
                if v.variant.is_frameshifting():
                    effective_variants.append(v)
        else:
            # If there are multiple frameshifts but the overall frames shifted is 0,
            # then checks if the reference sequence skipped has any stop codon.
            if len(upstream_indels) > 1:
                first_start = upstream_indels[0].variant.location.start
                last_end = upstream_indels[-1].variant.location.end
                rf_index = cds_start % 3
                if self.effect_analyzer.has_any_stop_codon_between(first_start, last_end, rf_index):
                    effective_variants.extend(upstream_indels)
        effective_variants.sort(key=lambda v: v.location)
        return effective_variants

    def call_peptides_main(self, variants:seqvar.VariantRecordPool,
            denylist:Set[str], check_variants:bool, check_canonical:bool,
            selenocysteine_termination:bool=False, is_mrna_end_nf:bool=False):
        """ Call peptide main """
        variant_peptides = set()
        tx_model = self.tx_model
        tx_seq = self.tx_seq
        rule = self.cleavage_params.enzyme
        exception = self.cleavage_params.exception
        chrom = tx_model.transcript.chrom
        codon_table = self.reference_data.codon_tables[chrom].codon_table
        start_codons = self.reference_data.codon_tables[chrom].start_codons

        is_circ_rna = False
        is_fusion = False

        is_coding = tx_model.is_protein_coding

        if variants[self.tx_id].fusion:
            is_fusion = True
            seq, variant_coordinates = self.sequence_builder.get_variant_sequence_fusion(
                seq=tx_seq.seq, variants=variants
            )
        elif variants[self.tx_id].circ_rna:
            is_circ_rna = True
            gene_seq = self.sequence_builder.get_gene_seq()
            seq, variant_coordinates = self.sequence_builder.get_variant_sequence_circ_rna(
                seq=gene_seq.seq, variants=variants
            )
        else:
            location = FeatureLocation(start=0, end=len(tx_seq.seq))
            seq, variant_coordinates = self.sequence_builder.get_variant_sequence(
                seq=tx_seq.seq, location=location, offset=0,
                variants=variants[self.tx_id].transcriptional, pool=variants
            )

        sec_positions = [] if is_circ_rna else \
            self.sequence_builder.get_sec_positions(variant_coordinates)
        variant_effects = self.effect_analyzer.check_variant_effect(str(seq), variant_coordinates)
        stop_lost, stop_gain, silent_mutation = variant_effects

        if not (is_coding and is_mrna_end_nf):
            cur_cds_end = len(seq)

        if not is_coding or is_circ_rna:
            alt_seq = dna.DNASeqRecord(seq)
            orf_tracker =ORFTracker(
                orf_starts=alt_seq.find_all_start_codons(start_codons)
            )
            if is_fusion:
                fusion_var = None
                # Get the last fusion, which should be the actual fusion
                # no insertions of intronic regions.
                for variant in reversed(variant_coordinates):
                    if variant.variant.is_fusion():
                        fusion_var = variant
                        break
                orf_tracker.filter_orfs_before(fusion_var.location.start)
        else:
            cds_start = tx_seq.orf.start
            orf_tracker = ORFTracker(
                orf_starts=[cds_start]
            )

        for cds_start in orf_tracker.orf_starts:
            cur_cds_end = len(seq) - (len(seq) - cds_start) % 3

            aa_seq = seq[cds_start:cur_cds_end].translate(table=codon_table, to_stop=False)
            if aa_seq[0] != 'M':
                aa_seq = Seq('M') + aa_seq[1:]
            if not is_circ_rna:
                for sec_start in sec_positions:
                    if (sec_start - cds_start) % 3 == 0:
                        sec_start_aa = int((sec_start - cds_start) / 3)
                        aa_seq = aa_seq[:sec_start_aa] + 'U' + aa_seq[sec_start_aa+1:]
            stop_index = aa_seq.find('*')
            if stop_index > -1:
                aa_seq = aa_seq[:stop_index]
            aa_seq = aa.AminoAcidSeqRecord(seq=aa_seq)

            sites = aa_seq.find_all_enzymatic_cleave_sites_with_ranges(rule, exception)

            # Finding the next M, so peptides that starts from the next M should
            # be processed with the correct `cds_start`
            if (not is_coding or is_circ_rna):
                next_inframe_cds = orf_tracker.get_next_in_frame_orf(cds_start)
                if next_inframe_cds == -1:
                    next_m = 0
                else:
                    next_m = (next_inframe_cds - cds_start) // 3
            else:
                next_m = 0

            sites.insert(0, (0, (0,1)))
            sites.append((len(aa_seq), (len(aa_seq), len(aa_seq))))
            for j, site in enumerate(sites[:-1]):
                lhs, lrange = site
                if 0 < next_m < lhs:
                    break
                last_inframe_cds = orf_tracker.get_last_in_frame_orf(cds_start, lhs)
                if last_inframe_cds > -1:
                    actual_cds_start = last_inframe_cds
                else:
                    actual_cds_start = cds_start

                # prev_cds_start = self.find_prev_cds_start_same_frame(
                #     cds_start=actual_cds_start,
                #     cds_start_positions=cds_start_positions
                # )

                for k in range(j + 1, min([j + 3, len(sites) - 1]) + 1):
                    rhs, rrange = sites[k]
                    tx_lhs = cds_start + lhs * 3
                    tx_rhs = cds_start + rhs * 3
                    if is_mrna_end_nf and tx_rhs + 3 > len(seq):
                        continue
                    peptide = aa_seq.seq[lhs:rhs]
                    is_in_denylist = str(peptide) in denylist \
                        and (not lhs == 0 or str(peptide[1:]) in denylist)
                    if is_in_denylist:
                        continue
                    tx_lrange = (cds_start + lrange[0] * 3, cds_start + lrange[1] * 3)
                    if k == len(sites) - 1 \
                            and tx_rhs + 3 <= len(seq) \
                            and seq[tx_rhs:tx_rhs + 3].translate() == '*':
                        tx_rrange = (tx_rhs, tx_rhs + 3)
                    else:
                        tx_rrange = (cds_start + rrange[0] * 3, cds_start + rrange[1] * 3)
                    if check_variants:
                        effective_variants = self.get_effective_variants(
                            lhs=tx_lhs, rhs=tx_rhs, lrange=tx_lrange, rrange=tx_rrange,
                            cds_start=actual_cds_start, variants=variant_coordinates,
                            variants_stop_lost=stop_lost, variants_stop_gain=stop_gain,
                            variants_silent_mutation=silent_mutation
                        )

                        if not effective_variants:
                            continue
                        if is_fusion and \
                                not any(v.variant.is_fusion() for v in effective_variants):
                            continue
                    else:
                        effective_variants = []

                    peptide_seqs = self.translational_modification(
                        peptide, lhs, tx_lhs, effective_variants, denylist,
                        check_variants, check_canonical, selenocysteine_termination
                    )
                    for peptide_seq in peptide_seqs:
                        variant_peptides.add(peptide_seq)
        return variant_peptides

    def translational_modification(self, seq:Seq, lhs:int, tx_lhs:int,
            effective_variants:List[VariantRecordWithCoordinate],
            denylist:List[str], check_variants:bool, check_canonical:bool,
            selenocysteine_termination:bool
            ) -> Iterable[str]:
        """ Apply any modification that could happen during translation. """
        candidates = []
        is_start = lhs == 0 and seq.startswith('M')
        if effective_variants or not check_variants:
            candidates.append(seq)
            if is_start:
                candidates.append(seq[1:])

        if selenocysteine_termination:
            k = 0
            while k > -1:
                k = seq.find('U', k)
                if k == -1:
                    break
                if not check_variants or \
                        any(v.location.start < tx_lhs + k * 3 for v in effective_variants):
                    candidates.append(seq[:k])
                    if is_start:
                        candidates.append(seq[1:k])
                k += 1

        # W > F
        if self.w2f:
            for candidate in copy.copy(candidates):
                w2f_indices = []
                i = 0
                while i > -1:
                    i = candidate.find('W', start=i)
                    if i > -1:
                        w2f_indices.append(i)
                        i += 1
                        if i > len(candidate):
                            break

                for k in range(1, len(w2f_indices) + 1):
                    for comb in combinations(w2f_indices, k):
                        seq_mod = candidate
                        for i in comb:
                            seq_mod_new = seq_mod[:i] + 'F'
                            if i + 1 < len(candidate):
                                seq_mod_new += seq_mod[i+1:]
                            seq_mod = seq_mod_new
                        candidates.append(seq_mod)

        for candidate in candidates:
            if self.peptide_is_valid(candidate, denylist, check_canonical):
                yield str(candidate)

    def generate_variant_comb(self, fusion:bool, circ_rna:bool
            ) -> Iterable[seqvar.VariantRecordPool]:
        """ Generate combination of variants. """
        variant_type_mapper:Dict[str, Tuple[seqvar.VariantRecord, str]] = {}
        start_index = self.tx_seq.orf.start + 3 if bool(self.tx_seq.orf) else 3
        for variant in self.variant_pool[self.tx_id].transcriptional:
            if variant.location.start == start_index - 1 \
                    and (variant.is_insertion() or variant.is_deletion()) \
                    and not variant.is_alternative_splicing():
                variant.to_end_inclusion(self.tx_seq)
            var_id = variant.get_minimal_identifier()
            variant_type_mapper[var_id] = (variant, 'transcriptional')
        for variant in self.variant_pool[self.tx_id].intronic:
            var_id = variant.get_minimal_identifier()
            variant_type_mapper[var_id] = (variant, 'intronic')
        if fusion:
            for variant in self.variant_pool[self.tx_id].fusion:
                if variant.location.start < start_index - 1:
                    continue
                var_id = variant.get_minimal_identifier()
                variant_type_mapper[var_id] = (variant, 'fusion')
                accepter_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
                if accepter_tx_id not in self.variant_pool:
                    continue
                accepter_var_series = self.variant_pool[accepter_tx_id]
                for accepter_var in accepter_var_series.transcriptional:
                    var_id = accepter_var.get_minimal_identifier()
                    variant_type_mapper[var_id] = (accepter_var, 'transcriptional')
                for accepter_var in accepter_var_series.intronic:
                    var_id = accepter_var.get_minimal_identifier()
                    variant_type_mapper[var_id] = (accepter_var, 'intronic')
        if circ_rna:
            for variant in self.variant_pool[self.tx_id].circ_rna:
                var_id = variant.get_minimal_identifier()
                variant_type_mapper[var_id] = (variant, 'circ_rna')

        all_variants = [v[0] for v in variant_type_mapper.values()]

        for i in range(len(all_variants)):
            for inds in combinations(range(len(all_variants)), i + 1):
                variants = [all_variants[i] for i in inds]
                if fusion \
                        and not any(variant_type_mapper[v.get_minimal_identifier()][1]
                                    == 'fusion' for v in variants):
                    continue
                if circ_rna \
                        and not any(variant_type_mapper[v.get_minimal_identifier()][1]
                                    == 'circ_rna' for v in variants):
                    continue
                pool = seqvar.VariantRecordPool()
                pool.anno = self.variant_pool.anno
                for variant in variants:
                    var_id = variant.get_minimal_identifier()
                    var_type = variant_type_mapper[var_id][1]
                    tx_id = variant.transcript_id
                    if var_type == 'transcriptional':
                        pool.add_transcriptional_variant(variant, tx_id)
                    elif var_type == 'intronic':
                        pool.add_intronic_variant(variant, tx_id)
                    elif var_type == 'fusion':
                        pool.add_fusion_variant(variant, tx_id)
                    elif var_type == 'circ_rna':
                        pool.add_circ_rna(variant, tx_id)
                pool.sort()
                if self.validator.has_incompatible_variants(pool):
                    continue
                yield pool


    def call_peptides(self):
        """ Call variant peptides """
        empty_pool = seqvar.VariantRecordPool()
        empty_pool[self.tx_id] = seqvar.TranscriptionalVariantSeries()
        denylist = self.call_peptides_main(
            variants=empty_pool, denylist=set(),
            check_variants=False, check_canonical=False,
            selenocysteine_termination=self.selenocysteine_termination,
            is_mrna_end_nf=False
        )
        main_peptides = set()

        for comb in self.generate_variant_comb(fusion=False, circ_rna=False):
            peptides = self.call_peptides_main(
                variants=comb, denylist=denylist,
                check_variants=True,  check_canonical=True,
                selenocysteine_termination=self.selenocysteine_termination,
                is_mrna_end_nf=self.tx_model.is_mrna_end_nf()
            )
            self.variant_peptides.update(peptides)
            main_peptides.update(peptides)

        for comb in self.generate_variant_comb(fusion=True, circ_rna=False):
            donor_tx_id = comb[self.tx_id].fusion[0].accepter_transcript_id
            donor_tx_model = self.reference_data.anno.transcripts[donor_tx_id]
            peptides = self.call_peptides_main(
                variants=comb, denylist=denylist,
                check_variants=True, check_canonical=True,
                is_mrna_end_nf=donor_tx_model.is_mrna_end_nf()
            )
            self.variant_peptides.update(peptides)

        denylist.update(main_peptides)
        for comb in self.generate_variant_comb(fusion=False, circ_rna=True):
            peptides = self.call_peptides_main(
                variants=comb, denylist=denylist,
                check_variants=True, check_canonical=True,
                is_mrna_end_nf=True
            )
            self.variant_peptides.update(peptides)

    def generate_variant_comb(self, fusion:bool, circ_rna:bool
            ) -> Iterable[seqvar.VariantRecordPool]:
        """ Generate combination of variants. """
        variant_type_mapper:Dict[str, Tuple[seqvar.VariantRecord, str]] = {}
        start_index = self.tx_seq.orf.start + 3 if bool(self.tx_seq.orf) else 3
        for variant in self.variant_pool[self.tx_id].transcriptional:
            if variant.location.start == start_index - 1 \
                    and (variant.is_insertion() or variant.is_deletion()) \
                    and not variant.is_alternative_splicing():
                variant.to_end_inclusion(self.tx_seq)
            var_id = variant.get_minimal_identifier()
            variant_type_mapper[var_id] = (variant, 'transcriptional')
        for variant in self.variant_pool[self.tx_id].intronic:
            var_id = variant.get_minimal_identifier()
            variant_type_mapper[var_id] = (variant, 'intronic')
        if fusion:
            for variant in self.variant_pool[self.tx_id].fusion:
                if variant.location.start < start_index - 1:
                    continue
                var_id = variant.get_minimal_identifier()
                variant_type_mapper[var_id] = (variant, 'fusion')
                accepter_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
                if accepter_tx_id not in self.variant_pool:
                    continue
                accepter_var_series = self.variant_pool[accepter_tx_id]
                for accepter_var in accepter_var_series.transcriptional:
                    var_id = accepter_var.get_minimal_identifier()
                    variant_type_mapper[var_id] = (accepter_var, 'transcriptional')
                for accepter_var in accepter_var_series.intronic:
                    var_id = accepter_var.get_minimal_identifier()
                    variant_type_mapper[var_id] = (accepter_var, 'intronic')
        if circ_rna:
            for variant in self.variant_pool[self.tx_id].circ_rna:
                var_id = variant.get_minimal_identifier()
                variant_type_mapper[var_id] = (variant, 'circ_rna')

        all_variants = [v[0] for v in variant_type_mapper.values()]

        for i in range(len(all_variants)):
            for inds in combinations(range(len(all_variants)), i + 1):
                variants = [all_variants[i] for i in inds]
                if fusion \
                        and not any(variant_type_mapper[v.get_minimal_identifier()][1]
                                    == 'fusion' for v in variants):
                    continue
                if circ_rna \
                        and not any(variant_type_mapper[v.get_minimal_identifier()][1]
                                    == 'circ_rna' for v in variants):
                    continue
                pool = seqvar.VariantRecordPool()
                pool.anno = self.variant_pool.anno
                for variant in variants:
                    var_id = variant.get_minimal_identifier()
                    var_type = variant_type_mapper[var_id][1]
                    tx_id = variant.transcript_id
                    if var_type == 'transcriptional':
                        pool.add_transcriptional_variant(variant, tx_id)
                    elif var_type == 'intronic':
                        pool.add_intronic_variant(variant, tx_id)
                    elif var_type == 'fusion':
                        pool.add_fusion_variant(variant, tx_id)
                    elif var_type == 'circ_rna':
                        pool.add_circ_rna(variant, tx_id)
                pool.sort()
                if self.validator.has_incompatible_variants(pool):
                    continue
                yield pool


