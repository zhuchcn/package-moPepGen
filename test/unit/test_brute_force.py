""" Unit tests for brute_force peptide window generation. """
import unittest
from Bio.Seq import Seq

from moPepGen import aa, params
from moPepGen.util.brute_force.peptide_site_generator import PeptideSiteGenerator


def _make_generator(miscleavage: int, min_len: int = 7, max_len: int = 25) -> PeptideSiteGenerator:
    cp = params.CleavageParams(
        enzyme='trypsin',
        exception=None,
        miscleavage=miscleavage,
        min_length=min_len,
        max_length=max_len,
        peptide_finding_mode='misc'
    )
    return PeptideSiteGenerator(cp)


class TestBruteForcePeptideSiteGenerator(unittest.TestCase):
    """Unit tests for PeptideSiteGenerator (brute_force)."""

    def test_misc_windows_miscleavage_combinations(self):
        """Misc mode should combine adjacent cleavage sites up to miscleavage.

        Sequence: MKAQRAA
        Trypsin cleavage sites at AA indices: after K (2) and after R (5)
        With sentinels 0 and 7, and miscleavage=1, expect windows:
          (0,2), (2,5), (5,7), (0,5), (2,7)
        """
        seq = aa.AminoAcidSeqRecord(seq=Seq('MKAQRAA'))
        gen = _make_generator(miscleavage=1)

        windows = gen.generate_peptide_windows(seq)
        pairs = {(lhs, rhs) for (lhs, rhs, _lrange, _rrange) in windows}

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

        windows = gen.generate_peptide_windows(seq)
        pairs = {(lhs, rhs) for (lhs, rhs, _lrange, _rrange) in windows}
        self.assertIn((0, 26), pairs)

