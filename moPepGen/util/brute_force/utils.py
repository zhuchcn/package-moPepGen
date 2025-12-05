""" Utility functions for brute force variant peptide caller """
from moPepGen import seqvar, params
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool


def create_mnvs(pool:VariantRecordPool, max_adjacent_as_mnv:int
        ) -> VariantRecordPool:
    """ Create MNVs """
    for tx_id in pool.data.keys():
        mnvs = seqvar.find_mnvs_from_adjacent_variants(
            pool[tx_id].transcriptional,
            max_adjacent_as_mnv
        )
        pool[tx_id].transcriptional = sorted(pool[tx_id].transcriptional + mnvs)

        mnvs = seqvar.find_mnvs_from_adjacent_variants(
            pool[tx_id].intronic,
            max_adjacent_as_mnv
        )
        pool[tx_id].intronic = sorted(pool[tx_id].intronic + mnvs)
    return pool

def fix_indel_after_start_codon(pool:VariantRecordPool,
        ref:params.ReferenceData) -> VariantRecordPool:
    """ Fix indel variants that are right after the start codon by shifting
    it to end inclusion to be consistant with callVariant """
    for tx_id in pool.data.keys():
        tx_model = ref.anno.transcripts[tx_id]
        chrom = tx_model.transcript.chrom
        tx_seq = tx_model.get_transcript_sequence(ref.genome[chrom])

        if tx_seq.orf:
            start = tx_seq.orf.start
        else:
            start = 0
        for v in pool[tx_id].transcriptional:
            if v.location.start == start + 2 \
                    and (v.is_insertion() or v.is_deletion()) \
                    and not v.is_fusion() \
                    and not v.is_alternative_splicing():
                v.to_end_inclusion(tx_seq)
    return pool
