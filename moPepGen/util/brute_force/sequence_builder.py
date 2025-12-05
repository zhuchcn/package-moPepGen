""" Sequence builder for variant sequences """
import copy
from typing import List, Tuple
from Bio.Seq import Seq
from moPepGen import seqvar, dna, constant
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool
from moPepGen.seqvar.VariantRecordWithCoordinate import VariantRecordWithCoordinate


class VariantSequenceBuilder:
    """
    Handles all variant sequence generation logic.
    
    Responsibilities:
    - Apply variants to sequences
    - Handle fusion variants
    - Handle circRNA variants
    - Track selenocysteine positions
    """
    
    def __init__(self, reference_data, tx_model, tx_seq, variant_pool, tx_id):
        self.reference_data = reference_data
        self.tx_model = tx_model
        self.tx_seq = tx_seq
        self.variant_pool = variant_pool
        self.tx_id = tx_id
        self.gene_seq = None  # Cached
    
    def get_gene_seq(self) -> dna.DNASeqRecord:
        """ Get the gene sequence and cache it if not already cached. """
        if self.gene_seq:
            return self.gene_seq
        gene_id = self.tx_model.gene_id
        gene_model = self.reference_data.anno.genes[gene_id]
        chrom = gene_model.chrom
        self.gene_seq = gene_model.get_gene_sequence(self.reference_data.genome[chrom])
        return self.gene_seq
    
    def get_variant_ref_seq(self, variant:VariantRecord) -> Seq:
        """ Get the reference sequence of a variant """
        if variant.type in ['Deletion', 'Substitution']:
            return self.tx_seq.seq[variant.location.start:variant.location.end]
        return variant.ref
    
    def get_variant_sequence(self, seq:Seq, location:FeatureLocation,
            offset:int, variants:List[VariantRecord],
            pool:VariantRecordPool
            ) -> Tuple[Seq, List[VariantRecordWithCoordinate]]:
        """ Get variant sequence by applying variants. """
        var_seq = seq
        variant_coordinates = []
        local_offset = 0
        for variant in variants:
            if variant.is_alternative_splicing():
                if variant.type == 'Deletion':
                    start = variant.location.start + local_offset - location.start
                    end = variant.location.end + local_offset - location.start
                    loc = FeatureLocation(start=start, end=start + 1)
                    variant_coordinate = VariantRecordWithCoordinate(
                        variant=variant,
                        location=loc
                    )
                    variant_coordinates.append(variant_coordinate)
                    alt_seq = var_seq[start:start+1]
                    local_offset = local_offset + len(alt_seq) - len(variant.location)
                    var_seq = var_seq[:start] + alt_seq + var_seq[end:]
                
                elif variant.type == 'Insertion':
                    start = variant.location.start + local_offset - location.start
                    end = variant.location.end + local_offset - location.start
                    
                    gene_seq = self.get_gene_seq()
                    donor_start = variant.get_donor_start()
                    donor_end = variant.get_donor_end()
                    alt_seq = str(gene_seq.seq[donor_start:donor_end])
                    loc = FeatureLocation(start=donor_start, end=donor_end)
                    insert_variants = [x for x in pool[self.tx_id].intronic
                        if loc.is_superset(x.location)]
                    alt_seq, insert_variants = self.get_variant_sequence(
                        seq=alt_seq, location=loc, offset=start,
                        variants=insert_variants, pool=pool
                    )
                    
                    variant_coordinate = VariantRecordWithCoordinate(
                        variant=variant,
                        location=FeatureLocation(start=start, end=start + len(alt_seq) + 1)
                    )
                    variant_coordinates.append(variant_coordinate)
                    variant_coordinates += insert_variants
                    local_offset = local_offset + len(alt_seq) + 1 - len(variant.location)
                    var_seq = var_seq[:start+1] + alt_seq + var_seq[end:]
                
                elif variant.type == 'Substitution':
                    start = variant.location.start + local_offset - location.start
                    end = variant.location.end + local_offset - location.start
                    
                    gene_seq = self.get_gene_seq()
                    donor_start = variant.get_donor_start()
                    donor_end = variant.get_donor_end()
                    alt_seq = str(gene_seq.seq[donor_start:donor_end])
                    loc = FeatureLocation(start=donor_start, end=donor_end)
                    insert_variants = [x for x in pool[self.tx_id].intronic
                        if loc.is_superset(x.location)]
                    alt_seq, insert_variants = self.get_variant_sequence(
                        seq=alt_seq, location=loc, offset=start,
                        variants=insert_variants, pool=pool
                    )
                    
                    variant_coordinate = VariantRecordWithCoordinate(
                        variant=variant,
                        location=FeatureLocation(start=start, end=start + len(alt_seq))
                    )
                    variant_coordinates.append(variant_coordinate)
                    variant_coordinates += insert_variants
                    local_offset = local_offset + len(alt_seq) - len(variant.location)
                    var_seq = var_seq[:start] + alt_seq + var_seq[end:]
            
            else:
                start = variant.location.start + local_offset - location.start
                end = variant.location.end + local_offset - location.start
                loc = FeatureLocation(
                    start=start + offset, end=start + len(variant.alt) + offset
                )
                variant_coordinate = VariantRecordWithCoordinate(
                    variant=variant,
                    location=loc
                )
                
                variant_coordinates.append(variant_coordinate)
                local_offset = local_offset + len(variant.alt) - len(variant.ref)
                var_seq = var_seq[:start] + variant.alt + var_seq[end:]
        
        return var_seq, variant_coordinates
    
    def get_variant_sequence_circ_rna(self, seq:Seq, variants:VariantRecordPool
            ) -> Tuple[Seq, List[VariantRecordWithCoordinate]]:
        """ Get the variant sequence of circRNA """
        number_of_circ = len(variants[self.tx_id].circ_rna)
        if number_of_circ != 1:
            raise ValueError(
                f"Should have exactly 1 circRNA, but {number_of_circ} were found."
            )
        circ = variants[self.tx_id].circ_rna[0]
        var_seq = Seq('')
        vars_coord = []
        
        for fragment in circ.fragments:
            loc = FeatureLocation(
                start=int(fragment.location.start), end=int(fragment.location.end)
            )
            frag_loc = FeatureLocation(
                start=int(fragment.location.start) + 3, end=int(fragment.location.end)
            )
            fragment = SeqFeature(location=frag_loc, chrom=fragment.chrom,
                attributes=fragment.attributes)
            new_seq = seq[loc.start:loc.end]
            frag_vars = variants.filter_variants(
                tx_ids=[circ.transcript_id], exclude_type=constant.ALTERNATIVE_SPLICING_TYPES,
                intron=False, segments=[fragment]
            )
            frag_seq, frag_vars_coord = self.get_variant_sequence(
                seq=new_seq, location=loc, offset=len(var_seq),
                variants=frag_vars, pool=variants
            )
            var_seq += frag_seq
            vars_coord += frag_vars_coord
        
        location = FeatureLocation(
            seqname=circ.gene_id,
            start=min(x.location.start for x in circ.fragments),
            end=max(x.location.end for x in circ.fragments)
        )
        circ_var = VariantRecord(
            location=location,
            ref=var_seq[0],
            alt='<circRNA>',
            _type='circRNA',
            _id=circ.id
        )
        
        vars_aloop = copy.deepcopy(vars_coord)
        seq_aloop = copy.deepcopy(var_seq)
        
        for _ in range(3):
            offset = len(var_seq)
            vars_extend = copy.deepcopy(vars_aloop)
            for variant in vars_extend:
                location = FeatureLocation(
                    seqname=variant.location.seqname,
                    start=variant.location.start + offset,
                    end=variant.location.end + offset
                )
                variant.location = location
            var_seq += seq_aloop
            vars_coord += vars_extend
        
        circ_var_coord = VariantRecordWithCoordinate(
            variant=circ_var, location=FeatureLocation(start=0, end=len(var_seq))
        )
        vars_coord.insert(0, circ_var_coord)
        
        return var_seq, vars_coord
    
    def get_variant_sequence_fusion(self, seq:Seq, variants:VariantRecordPool
            ) -> Tuple[Seq, List[VariantRecordWithCoordinate]]:
        """ Get variant sequence with fusion. """
        number_of_fusion = len(variants[self.tx_id].fusion)
        if not number_of_fusion == 1:
            raise ValueError(
                f"Should have exactly 1 fusion, but {number_of_fusion} were found."
            )
        fusion = variants[self.tx_id].fusion[0]
        var_seq = seq[:fusion.location.start]
        location = FeatureLocation(start=0, end=len(var_seq))
        var_seq, variant_coordinates = self.get_variant_sequence(
            seq=var_seq, location=location, offset=0,
            variants=variants[self.tx_id].transcriptional, pool=variants
        )
        
        left_insert_start = fusion.attrs['LEFT_INSERTION_START']
        left_insert_end = fusion.attrs['LEFT_INSERTION_END']
        right_insert_start = fusion.attrs['RIGHT_INSERTION_START']
        right_insert_end = fusion.attrs['RIGHT_INSERTION_END']
        right_tx_id = fusion.attrs['ACCEPTER_TRANSCRIPT_ID']
        right_gene_id = fusion.attrs['ACCEPTER_GENE_ID']
        
        additional_seq = Seq('')
        additional_variants:List[VariantRecordWithCoordinate] = []
        
        # left insertion
        if left_insert_start is not None:
            gene_seq = self.get_gene_seq()
            location = FeatureLocation(start=left_insert_start, end=left_insert_end)
            insert_seq = gene_seq.seq[left_insert_start:left_insert_end]
            insert_variants = [x for x in variants[self.tx_id].intronic
                if location.is_superset(x.location)]
            insert_seq, insert_variants = self.get_variant_sequence(
                seq=insert_seq, location=location,
                offset=len(var_seq) + len(additional_seq),
                variants=insert_variants, pool=variants
            )
            location=FeatureLocation(
                start=len(var_seq),
                end=len(var_seq) + len(insert_seq)
            )
            insert_fusion = VariantRecordWithCoordinate(
                variant=fusion, location=location
            )
            insert_variants.insert(0, insert_fusion)
            additional_seq += insert_seq
            additional_variants += insert_variants
        
        # right insertion
        if right_insert_start is not None:
            gene_model = self.reference_data.anno.genes[right_gene_id]
            chrom = gene_model.chrom
            gene_seq = gene_model.get_gene_sequence(self.reference_data.genome[chrom])
            location = FeatureLocation(start=right_insert_start, end=right_insert_end)
            insert_seq = gene_seq.seq[right_insert_start:right_insert_end]
            insert_seq, insert_variants = self.get_variant_sequence(
                seq=insert_seq, location=location,
                offset=len(var_seq) + len(additional_seq),
                variants=variants[right_tx_id].intronic if right_tx_id in variants else [],
                pool=variants
            )
            location=FeatureLocation(
                start=len(var_seq) + len(additional_seq),
                end=len(var_seq) + len(additional_seq) + len(insert_seq)
            )
            insert_fusion = VariantRecordWithCoordinate(
                variant=fusion, location=location
            )
            insert_variants.insert(0, insert_fusion)
            additional_seq += insert_seq
            additional_variants += insert_variants
        
        right_tx_model = self.reference_data.anno.transcripts[right_tx_id]
        accepter_chrom = right_tx_model.transcript.location.seqname
        breakpoint_gene = fusion.get_accepter_position()
        breakpoint_tx = self.reference_data.anno.coordinate_gene_to_transcript(
            index=breakpoint_gene,
            gene=right_gene_id,
            transcript=right_tx_id
        )
        accepter_seq = right_tx_model.get_transcript_sequence(
            self.reference_data.genome[accepter_chrom]
        )
        location = FeatureLocation(start=breakpoint_tx, end=len(accepter_seq))
        insert_seq = accepter_seq.seq[breakpoint_tx:]
        insert_seq, insert_variants = self.get_variant_sequence(
            seq=insert_seq, location=location,
            offset=len(var_seq) + len(additional_seq),
            variants=variants[right_tx_id].transcriptional if right_tx_id in variants else [],
            pool=variants
        )
        
        location = FeatureLocation(
            start=len(var_seq) + len(additional_seq),
            end=len(var_seq) + len(additional_seq) + len(insert_seq)
        )
        fusion_var = VariantRecordWithCoordinate(variant=fusion, location=location)
        insert_variants.insert(0, fusion_var)
        
        additional_seq += insert_seq
        additional_variants += insert_variants
        
        var_seq += additional_seq
        variant_coordinates += additional_variants
        
        return var_seq, variant_coordinates
    
    def get_sec_positions(self, variants:List[VariantRecordWithCoordinate]) -> List[int]:
        """ Get Sec positions in the altered sequence. """
        sec_iter = iter(self.tx_seq.selenocysteine)
        var_iter = iter(variants)
        sec_i = next(sec_iter, None)
        var_i = next(var_iter, None)
        
        sec_positions = []
        offset = 0
        fusion_breakpoint = None
        while sec_i:
            if var_i:
                # Not consider Sec after fusion breakpoint.
                if var_i.variant.is_fusion():
                    fusion_breakpoint = var_i.variant.location.start
                    var_i = None
                    continue
                if var_i.variant.location.end <= sec_i.start:
                    ref_len = var_i.variant.get_ref_len()
                    alt_len = var_i.variant.get_alt_len()
                    offset += (alt_len - ref_len)
                    var_i = next(var_iter, None)
                    continue
                
                if var_i.variant.type == 'Deletion':
                    # Because for Deletion, the sequence after the REF
                    # nucleotide is deleted so the first nucleotide is unchanged.
                    var_loc = FeatureLocation(
                        var_i.variant.location.start+1,
                        var_i.variant.location.end
                    )
                else:
                    var_loc = var_i.variant.location
                
                if var_loc.overlaps(sec_i):
                    sec_i = next(sec_iter, None)
                    continue
            if fusion_breakpoint and sec_i.end >= fusion_breakpoint:
                break
            sec_start = sec_i.start + offset
            sec_positions.append(sec_start)
            sec_i = next(sec_iter, None)
        return sec_positions
