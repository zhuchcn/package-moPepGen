""" Peptide window generator for peptide finding modes. """
from __future__ import annotations
from typing import TYPE_CHECKING
from moPepGen import constant

if TYPE_CHECKING:
    from typing import List, Tuple
    from moPepGen import aa, params


class PeptideSiteGenerator:
    """Generate concrete peptide windows for each peptide-finding mode.

    The caller no longer loops over cleavage *sites* to form windows. Instead
    this class returns explicit (lhs, rhs, lrange, rrange) tuples that already
    honor miscleavage rules (for misc) or the sliding/archipel window lengths.
    """

    def __init__(self, cleavage_params: 'params.CleavageParams'):
        self.mode = cleavage_params.peptide_finding_mode
        self.enzyme = cleavage_params.enzyme
        self.exception = cleavage_params.exception
        self.flanking_size = cleavage_params.flanking_size
        self.miscleavage = cleavage_params.miscleavage
        self.min_length = cleavage_params.min_length
        self.max_length = cleavage_params.max_length

    def generate_peptide_windows(self, aa_seq: 'aa.AminoAcidSeqRecord'
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Return peptide windows as (lhs, rhs, lrange, rrange).

        - misc: Uses enzymatic cleavage sites and combines left/right sites
          allowing up to `miscleavage` missed cleavages.
        - sliding_window: All start positions with lengths in [min_length, max_length].
        - archipel: Same enumeration as sliding_window; downstream filtering uses
          variant context and flanking size.
        """
        if self.mode == constant.PeptideFindingMode.MISC.value:
            return self._generate_misc_windows(aa_seq)
        if self.mode == constant.PeptideFindingMode.SLIDING_WINDOW.value:
            return self._generate_sliding_windows(aa_seq)
        if self.mode == constant.PeptideFindingMode.ARCHIPEL.value:
            return self._generate_archipel_windows(aa_seq)
        raise ValueError(f"Unknown peptide finding mode: {self.mode}")

    def _generate_misc_windows(self, aa_seq: 'aa.AminoAcidSeqRecord'
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Generate cleavage-based windows honoring miscleavage allowance."""
        sites = aa_seq.find_all_enzymatic_cleave_sites_with_ranges(
            self.enzyme, self.exception
        )

        # Ensure start/end sentinels are present
        sites = [(0, (0, 1))] + sites + [(len(aa_seq), (len(aa_seq), len(aa_seq)))]

        windows: List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]] = []
        max_sites = len(sites)
        for i, (lhs, lrange) in enumerate(sites[:-1]):
            max_j = min(i + self.miscleavage + 1, max_sites - 1)
            for j in range(i + 1, max_j + 1):
                rhs, rrange = sites[j]
                # Do not filter by length here: caller may trim leading M
                # (producing a shorter valid peptide). Final filters are applied
                # after translational modifications.
                windows.append((lhs, rhs, lrange, rrange))
        return windows

    def _generate_sliding_windows(self, aa_seq: 'aa.AminoAcidSeqRecord'
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Generate fixed-length windows for sliding window mode."""
        windows: List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]] = []
        max_len = min(int(self.max_length), len(aa_seq))
        min_len = self.min_length
        if min_len > max_len:
            return windows
        for lhs in range(0, len(aa_seq)):
            for length in range(min_len, max_len + 1):
                rhs = lhs + length
                if rhs > len(aa_seq):
                    break
                # Minimal ranges that align with how caller expects lrange/rrange
                windows.append((lhs, rhs, (lhs, lhs + 1), (rhs, rhs + 1)))
        return windows

    def _generate_archipel_windows(self, aa_seq: 'aa.AminoAcidSeqRecord'
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Generate windows for archipel mode (same enumeration as sliding)."""
        return self._generate_sliding_windows(aa_seq)
