""" Variant compatibility validator """
from typing import List, Dict
from moPepGen import seqvar, constant
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool
from moPepGen.seqvar.VariantRecordWithCoordinate import VariantRecordWithCoordinate


class VariantCompatibilityValidator:
    """
    Validates variant combinations for compatibility.

    Responsibilities:
    - Check for overlapping variants
    - Validate fusion/splicing variant compatibility
    - Check variant combinations
    - Validate circRNA variants
    """

    def __init__(self, reference_data, tx_model, tx_seq, variant_pool, tx_id):
        self.reference_data = reference_data
        self.tx_model = tx_model
        self.tx_seq = tx_seq
        self.variant_pool = variant_pool
        self.tx_id = tx_id

    @staticmethod
    def has_overlapping_variants(variants:List[VariantRecord]) -> bool:
        """ Checks if any variants overlap. """
        for i, left in enumerate(variants):
            if i == len(variants) - 1:
                continue
            for right in variants[i+1:]:
                if left.location.end >= right.location.start:
                    return True
        return False

    def should_clip_trailing_nodes(self, variants:List[VariantRecordWithCoordinate]):
        """ Checks whether the trailing nodes should be excluded """
        return any(v.variant.is_circ_rna() for v in variants) \
            or self.tx_model.is_mrna_end_nf()

    def has_any_invalid_variants_on_inserted_sequences(self,
            pool:VariantRecordPool) -> bool:
        """ Checks if any variants carried by fusion or alt splicing. Invalid
        variants are those not in the region of sequence introduced by fusion.
        For example, a variant of the donor transcript before the breakpoint is
        invalid.
        """
        alt_splices = [x for x in pool[self.tx_id].transcriptional if x.is_alternative_splicing()]
        if not pool[self.tx_id].fusion and not alt_splices:
            return False

        inserted_intronic_region:Dict[str, List[FeatureLocation]] = {}

        if pool[self.tx_id].fusion:
            fusion = pool[self.tx_id].fusion[0]
            left_insert_start = fusion.attrs['LEFT_INSERTION_START']
            left_insert_end = fusion.attrs['LEFT_INSERTION_END']
            right_insert_start = fusion.attrs['RIGHT_INSERTION_START']
            right_insert_end = fusion.attrs['RIGHT_INSERTION_END']
            right_tx_id = fusion.attrs['ACCEPTER_TRANSCRIPT_ID']
            right_tx_model = self.reference_data.anno.transcripts[right_tx_id]
            right_gene_id = right_tx_model.gene_id
            right_breakpoint_gene = fusion.get_accepter_position()
            right_breakpoint_tx = self.reference_data.anno.coordinate_gene_to_transcript(
                index=right_breakpoint_gene,
                gene=right_gene_id,
                transcript=right_tx_id
            )

            if right_tx_id in pool:
                alt_splices += [x for x in pool[right_tx_id].transcriptional
                    if x.is_alternative_splicing()]

            if any(x.location.end >= fusion.location.start
                    for x in pool[self.tx_id].transcriptional):
                return True

            if left_insert_start:
                loc = FeatureLocation(start=left_insert_start, end=left_insert_end)
                if self.tx_id not in inserted_intronic_region:
                    inserted_intronic_region[self.tx_id] = []
                inserted_intronic_region[self.tx_id].append(loc)

            if right_insert_start:
                loc = FeatureLocation(start=right_insert_start, end=right_insert_end)
                if right_tx_id not in inserted_intronic_region:
                    inserted_intronic_region[right_tx_id] = []
                inserted_intronic_region[right_tx_id].append(loc)

            if right_tx_id in pool \
                    and any(x.location.start <= right_breakpoint_tx
                        for x in pool[right_tx_id].transcriptional):
                return True

        for alt_splice in alt_splices:
            donor_start = alt_splice.attrs.get('DONOR_START')
            if not donor_start:
                continue
            donor_start = int(donor_start)
            donor_end = int(alt_splice.attrs['DONOR_END'])
            loc = FeatureLocation(start=donor_start, end=donor_end)
            if alt_splice.transcript_id not in inserted_intronic_region:
                inserted_intronic_region[alt_splice.transcript_id] = []
            inserted_intronic_region[alt_splice.transcript_id].append(loc)

        if not inserted_intronic_region:
            return False

        for tx_id in pool:
            if not pool[tx_id].intronic:
                continue
            if tx_id not in inserted_intronic_region:
                return True
            for v in pool[tx_id].intronic:
                if not any(x.start < v.location.start and x.end > v.location.end
                        for x in inserted_intronic_region[tx_id]):
                    return True
        return False

    @staticmethod
    def has_any_invalid_variant_on_circ() -> bool:
        """ Checks if any variants are invalid with circRNA. """
        return False

    def has_incompatible_variants(self, pool:VariantRecordPool) -> bool:
        """ Whether there are incompatible variants, i.e. variants that overlap
        with each other. """
        # check if there is any variants
        if self.tx_id not in pool:
            return True

        if not pool[self.tx_id].transcriptional\
                and not pool[self.tx_id].circ_rna \
                and not pool[self.tx_id].fusion:
            return True

        # check if there are multiple novel transcript variants (fusion + circ)
        n_fusion = len(pool[self.tx_id].fusion)
        n_circ = len(pool[self.tx_id].circ_rna)
        n_alt_splice = len([x for x in pool[self.tx_id].transcriptional
            if x.type in constant.ALTERNATIVE_SPLICING_TYPES])

        orf = self.tx_seq.orf
        is_mrna_end_nf = self.tx_model.is_mrna_end_nf()
        start_index = orf.start + 3 if bool(orf) else 3

        if n_circ == 0:
            for variant in pool[self.tx_id].transcriptional:
                if variant.location.start < start_index:
                    return True
                if is_mrna_end_nf and orf is not None:
                    orf_end_trinuc = FeatureLocation(start=orf.end-3, end=orf.end)
                    if variant.location.overlaps(orf_end_trinuc):
                        return True

        for fusion in pool[self.tx_id].fusion:
            if fusion.location.start < start_index:
                return True
            accepter_tx_id = fusion.attrs['ACCEPTER_TRANSCRIPT_ID']
            if accepter_tx_id not in pool:
                continue
            n_alt_splice += len([x for x in pool[accepter_tx_id].transcriptional
                if x.type in constant.ALTERNATIVE_SPLICING_TYPES])
        if n_fusion + n_circ > 1:
            return True

        # not allowing any alternative splicing comb with fusion or circ
        if (n_fusion > 0 or n_circ > 0) and n_alt_splice > 0:
            return True

        if self.has_any_invalid_variants_on_inserted_sequences(pool):
            return True

        # if self.has_any_invalid_variant_on_circ(pool):
        #     return True

        for series in pool.data.values():
            if self.has_overlapping_variants(series.transcriptional):
                return True
            if self.has_overlapping_variants(series.intronic):
                return True

        if self.tx_model.is_mrna_end_nf():
            for variant in pool[self.tx_id].transcriptional:
                if variant.location.end >= self.tx_seq.orf.end:
                    return True
                if n_circ == 0 \
                        and variant.location.start >= self.tx_seq.orf.end - 3:
                    return True
        return False
