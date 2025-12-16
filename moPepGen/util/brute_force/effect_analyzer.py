""" Variant effect analyzer """
import math
from typing import List, Tuple
from Bio.Seq import Seq
from moPepGen import seqvar
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.seqvar.VariantRecordWithCoordinate import VariantRecordWithCoordinate


class VariantEffectAnalyzer:
    """
    Analyzes variant effects on translation.

    Responsibilities:
    - Check stop loss
    - Check stop gain
    - Check silent mutations
    - Check for stop codons in range
    """

    def __init__(self, tx_seq, tx_model, tx_id):
        self.tx_seq = tx_seq
        self.tx_model = tx_model
        self.tx_id = tx_id

    def is_stop_lost(self, variant:VariantRecord, reading_frame_index:int) -> bool:
        """ Check whether the variant is a stop lost mutation. """
        orf_index = (self.tx_seq.orf.start if self.tx_seq.orf else 0) % 3

        if self.tx_model.is_protein_coding:
            if orf_index != reading_frame_index:
                return False
            orf_end = self.tx_seq.orf.end
            stop_codon = FeatureLocation(orf_end, orf_end + 3)
            return variant.location.overlaps(stop_codon)

        i = variant.variant.location.start - (variant.location.start - orf_index) % 3
        while i < variant.location.end:
            if self.tx_seq.seq[i:i+3] in ['TAA', 'TAG', 'TGA']:
                return True
            i += 3
        return False

    def has_any_stop_codon_between(self, start:int, end:int, rf_index:int) -> bool:
        """ Checks if there is any stop codon between a given range of the tx. """
        aa_start = math.floor((start - rf_index) / 3)
        aa_end = math.ceil((end - rf_index) / 3)
        aa_end = min(aa_end, math.floor((len(self.tx_seq) - rf_index)/3))
        aa_seq:Seq = self.tx_seq.seq[rf_index:].translate()[aa_start:aa_end]
        i = 0
        sec_sites = {int(x.start) for x in self.tx_seq.selenocysteine}
        while i > -1:
            i = aa_seq.find('*', i)
            if i == -1:
                return False
            if aa_start * 3 + rf_index + i * 3 not in sec_sites:
                return True
            i += 1
        return False

    def check_variant_effect(self, seq:str,
            variants:List[VariantRecordWithCoordinate]
            ) -> Tuple[List[Tuple[bool, bool, bool]],
                      List[Tuple[bool, bool, bool]],
                      List[Tuple[bool, bool, bool]]]:
        """ Check the variant effects, including stop lost, stop gain, and
        silent mutation """
        stop_lost:List[Tuple[bool, bool, bool]] = []
        stop_gain:List[Tuple[bool, bool, bool]] = []
        silent_mutation:List[Tuple[bool, bool, bool]] = []

        for variant in variants:
            skip_stop_lost = False
            skip_stop_gain = False
            skip_silent_mutation = False
            stop_gain_i:List[bool] = []
            stop_lost_i:List[bool] = []
            silent_mutation_i:List[bool] = []

            if variant.variant.is_fusion() or variant.variant.is_circ_rna():
                stop_gain_i = (False, False, False)
                skip_stop_gain = True

            if variant.variant.is_fusion() \
                    or variant.variant.is_circ_rna() \
                    or variant.variant.transcript_id != self.tx_id \
                    or 'TRANSCRIPT_ID' in variant.variant.attrs:
                stop_lost_i = (False, False, False)
                skip_stop_lost = True

            if variant.variant.is_fusion() \
                    or variant.variant.is_circ_rna() \
                    or variant.variant.is_alternative_splicing()\
                    or variant.variant.is_frameshifting():
                silent_mutation_i = (False, False, False)
                skip_silent_mutation = True

            for i in range(3):
                lhs_offset = (variant.location.start - i) % 3
                alt_lhs = variant.location.start - lhs_offset
                var_seq = seq[alt_lhs:variant.location.end]
                n_carry_over = (3 - (len(var_seq) % 3)) % 3
                alt_rhs = min(len(seq), variant.location.end + n_carry_over)
                var_seq += seq[variant.location.end:alt_rhs]

                loc = variant.variant.location
                lhs = loc.start - lhs_offset
                rhs = loc.end + (3 - (loc.end - lhs) % 3) % 3
                rhs = min(len(self.tx_seq), rhs)
                ref_seq = self.tx_seq.seq[lhs:rhs]
                ref_rf_index = lhs % 3

                var_aa = Seq(var_seq).translate(to_stop=False)
                ref_aa = Seq(ref_seq).translate(to_stop=False)
                for sec in self.tx_seq.selenocysteine:
                    if lhs <= sec.start < sec.end <= rhs \
                            and (sec.start - lhs) % 3 == 0:
                        sec_aa = int((sec.start - lhs) / 3)
                        ref_aa = ref_aa[:sec_aa] + 'U' + ref_aa[sec_aa+1:]

                    # this is when the variant is a deletion, the position is
                    # after the 3rd nt of selenocysteine codon, the codon sequence
                    # is unchanged and the amino acid sequence should be replaced
                    # to U from *.
                    if variant.variant.is_deletion():
                        if not variant.variant.is_end_inclusion() \
                                and lhs_offset == 2 \
                                and alt_lhs == sec.start and sec.end == alt_rhs \
                                and (sec.start - alt_lhs) % 3 == 0:
                            sec_aa = int((sec.start - alt_lhs) / 3)
                            var_aa = var_aa[:sec_aa] + 'U' + var_aa[sec_aa+1:]

                if not skip_stop_lost:
                    is_stop_lost = '*' not in var_aa and '*' in ref_aa \
                        and (self.tx_seq.orf is None
                            or ref_rf_index == self.tx_seq.orf.start % 3)
                    stop_lost_i.append(is_stop_lost)

                sec_altering = any(loc.overlaps(sec) for sec in self.tx_seq.selenocysteine)

                if not skip_stop_gain:
                    is_stop_gain = var_aa == '' or var_aa.startswith('*') \
                        and (
                            not (ref_aa == '' or ref_aa.startswith('*'))
                            or sec_altering
                        )
                    stop_gain_i.append(is_stop_gain)

                if not skip_silent_mutation:
                    is_silent_mutation = ref_aa == var_aa and not sec_altering
                    silent_mutation_i.append(is_silent_mutation)

            stop_lost.append(tuple(stop_lost_i))
            stop_gain.append(tuple(stop_gain_i))
            silent_mutation.append(tuple(silent_mutation_i))

        return stop_lost, stop_gain, silent_mutation
