""" Unit tests for brute_force peptide window generation. """
import unittest
from Bio.Seq import Seq

from moPepGen import aa, params
from moPepGen.util.brute_force.peptide_candidate_generator import PeptideCandidateGenerator


class _DummyTracker:
    """ Dummy ORFTracker for testing. """
    def get_next_in_frame_orf(self, _cds_start:int):
        """ Get the next in-frame ORF. """
        return -1

    def get_last_in_frame_orf(self, _cds_start:int, _lhs:int):
        """ Get the last in-frame ORF before lhs. """
        return -1

class _DummyEffectAnalyzer:
    """ Dummy VariantEffectAnalyzer for testing. """
    def has_any_stop_codon_between(self, *_args, **_kwargs):
        """ Check if there is any stop codon between positions. """
        return False


def _make_generator(miscleavage: int, min_len: int = 7, max_len: int = 25
        ) -> PeptideCandidateGenerator:
    """ Create a minimal PeptideCandidateGenerator for misc mode testing. """
    cp = params.CleavageParams(
        enzyme='trypsin',
        exception=None,
        miscleavage=miscleavage,
        min_length=min_len,
        max_length=max_len,
        peptide_finding_mode='misc'
    )
    # Minimal generator for misc mode tests: no variants, not coding, no fusion/circ.
    return PeptideCandidateGenerator(
        cleavage_params=cp,
        seq=Seq(''),
        is_coding=False,
        is_fusion=False,
        is_circ_rna=False,
        is_mrna_end_nf=False,
        check_variants=False,
        variants=[],
        stop_lost=[],
        stop_gain=[],
        silent_mutation=[],
        denylist=set(),
        orf_tracker=_DummyTracker(),
        effect_analyzer=_DummyEffectAnalyzer(),
    )


class TestBruteForcePeptideCandidateGenerator(unittest.TestCase):
    """Unit tests for PeptideCandidateGenerator (brute_force)."""

    def test_misc_windows_miscleavage_combinations(self):
        """Misc mode should combine adjacent cleavage sites up to miscleavage.

        Sequence: MKAQRAA
        Trypsin cleavage sites at AA indices: after K (2) and after R (5)
        With sentinels 0 and 7, and miscleavage=1, expect windows:
          (0,2), (2,5), (5,7), (0,5), (2,7)
        """
        seq = aa.AminoAcidSeqRecord(seq=Seq('MKAQRAA'))
        gen = _make_generator(miscleavage=1)

        candidates = list(gen.generate_peptide_candidates(seq, cds_start=0))
        pairs = {(c.lhs, c.rhs) for c in candidates}

        expected = {(0, 2), (2, 5), (5, 7), (0, 5), (2, 7)}
        self.assertEqual(pairs, expected)

    def test_misc_windows_include_overlength_for_startM_trimming(self):
        """Misc mode should include overlength window to allow start-M trimming.

        Sequence of length 26 with no tryptic sites → only window (0, 26).
        Even if max_length=25, the (0,26) window must be present so that the
        caller can trim the initial M and then enforce length constraints later.
        """
        seq = aa.AminoAcidSeqRecord(seq=Seq('M' + 'A' * 25))
        gen = _make_generator(miscleavage=2, min_len=7, max_len=25)

        candidates = list(gen.generate_peptide_candidates(seq, cds_start=0))
        pairs = {(c.lhs, c.rhs) for c in candidates}
        self.assertIn((0, 26), pairs)
