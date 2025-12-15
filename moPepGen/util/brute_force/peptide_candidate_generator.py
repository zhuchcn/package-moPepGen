""" Peptide window generator for peptide finding modes. """
from __future__ import annotations
from typing import TYPE_CHECKING
import dataclasses
from moPepGen import constant
from moPepGen.SeqFeature import FeatureLocation

if TYPE_CHECKING:
    from typing import List, Tuple, Set, Iterable
    from Bio.Seq import Seq
    from moPepGen.aa import AminoAcidSeqRecord
    from moPepGen.seqvar.VariantRecordWithCoordinate import VariantRecordWithCoordinate
    from moPepGen.params import CleavageParams
    from .orf_tracker import ORFTracker
    from .effect_analyzer import VariantEffectAnalyzer


@dataclasses.dataclass
class PeptideCandidate:
    """Data class for a peptide candidate window."""
    peptide: str
    lhs: int
    rhs: int
    lhs_tx: int
    rhs_tx: int
    lrange: Tuple[int, int]
    rrange: Tuple[int, int]
    effective_variants: List[VariantRecordWithCoordinate]

class PeptideCandidateGenerator:
    """Generate concrete peptide candidates for each peptide-finding mode.

    This class generates explicit (lhs, rhs, lrange, rrange) window tuples that
    honor miscleavage rules (for misc), sliding/archipel length constraints,
    and variant context (for archipel).
    """

    def __init__(self, cleavage_params:CleavageParams, seq:Seq, is_coding:bool,
            is_fusion:bool, is_circ_rna:bool, is_mrna_end_nf:bool, check_variants:bool,
            variants:List[VariantRecordWithCoordinate],
            stop_lost:List[Tuple[bool,bool,bool]],
            stop_gain:List[Tuple[bool,bool,bool]],
            silent_mutation:List[Tuple[bool,bool,bool]],
            denylist:Set[str], orf_tracker:ORFTracker,
            effect_analyzer:VariantEffectAnalyzer):
        self.cleavage_params = cleavage_params
        self.seq = seq
        self.is_coding = is_coding
        self.is_fusion = is_fusion
        self.is_circ_rna = is_circ_rna
        self.is_mrna_end_nf = is_mrna_end_nf
        self.denylist = denylist
        self.check_variants = check_variants
        self.variants = variants
        self.stop_lost = stop_lost
        self.stop_gain = stop_gain
        self.silent_mutation = silent_mutation
        self.orf_tracker = orf_tracker
        self.effect_analyzer = effect_analyzer

    def generate_peptide_candidates(self, aa_seq: AminoAcidSeqRecord,
            cds_start: int = None
            ) -> Iterable[PeptideCandidate]:
        """Return peptide candidates as (lhs, rhs, lrange, rrange).

        - misc: Uses enzymatic cleavage sites and combines left/right sites
          allowing up to `miscleavage` missed cleavages.
        - sliding_window: All start positions with lengths in [min_length, max_length].
        - archipel: Same enumeration as sliding_window; downstream filtering uses
          variant context and flanking size.
        """
        mode = self.cleavage_params.peptide_finding_mode
        if mode == constant.PeptideFindingMode.MISC.value:
            windows = self._generate_misc_windows(aa_seq)
        elif mode == constant.PeptideFindingMode.SLIDING_WINDOW.value:
            windows = self._generate_sliding_windows(aa_seq)
        elif mode == constant.PeptideFindingMode.ARCHIPEL.value:
            windows =  self._generate_archipel_windows(aa_seq, cds_start)
        else:
            raise ValueError(f"Unknown peptide finding mode: {mode}")

        # Finding the next M, so peptides that starts from the next M should
        # be processed with the correct `cds_start`
        if (not self.is_coding or self.is_circ_rna):
            next_inframe_cds = self.orf_tracker.get_next_in_frame_orf(cds_start)
            if next_inframe_cds == -1:
                next_m = 0
            else:
                next_m = (next_inframe_cds - cds_start) // 3
        else:
            next_m = 0

        for lhs, rhs, lrange, rrange in windows:
            if 0 < next_m < lhs:
                break
            # Find the appropriate CDS start for variant effect evaluation.
            # When multiple in-frame starts exist, we need to anchor variant
            # checks to the most recent upstream start, not the current loop's
            # cds_start.
            last_inframe_cds = self.orf_tracker.get_last_in_frame_orf(cds_start, lhs)
            actual_cds_start = last_inframe_cds if last_inframe_cds > -1 else cds_start

            tx_lhs = cds_start + lhs * 3
            tx_rhs = cds_start + rhs * 3
            if self.is_mrna_end_nf and tx_rhs + 3 > len(self.seq):
                continue
            peptide = aa_seq.seq[lhs:rhs]
            is_in_denylist = str(peptide) in self.denylist \
                and (lhs != 0 or str(peptide[1:]) in self.denylist)
            if is_in_denylist:
                continue
            tx_lrange = (cds_start + lrange[0] * 3, cds_start + lrange[1] * 3)
            if rhs == len(aa_seq) \
                    and tx_rhs + 3 <= len(self.seq) \
                    and self.seq[tx_rhs:tx_rhs + 3].translate() == '*':
                tx_rrange = (tx_rhs, tx_rhs + 3)
            else:
                tx_rrange = (cds_start + rrange[0] * 3, cds_start + rrange[1] * 3)
            if self.check_variants:
                effective_variants = self.get_effective_variants(
                    lhs=tx_lhs, rhs=tx_rhs, lrange=tx_lrange, rrange=tx_rrange,
                    cds_start=actual_cds_start, variants=self.variants,
                    variants_stop_lost=self.stop_lost, variants_stop_gain=self.stop_gain,
                    variants_silent_mutation=self.silent_mutation
                )

                if not effective_variants:
                    continue
                if self.is_fusion and \
                        not any(v.variant.is_fusion() for v in effective_variants):
                    continue
            else:
                effective_variants = []

            yield PeptideCandidate(
                peptide=peptide,
                lhs=lhs,
                rhs=rhs,
                lhs_tx=tx_lhs,
                rhs_tx=tx_rhs,
                lrange=lrange,
                rrange=rrange,
                effective_variants=effective_variants
            )

    def _generate_misc_windows(self, aa_seq: AminoAcidSeqRecord
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Generate cleavage-based windows honoring miscleavage allowance."""
        sites = aa_seq.find_all_enzymatic_cleave_sites_with_ranges(
            rule=self.cleavage_params.enzyme,
            exception=self.cleavage_params.exception
        )

        # Ensure start/end sentinels are present
        sites = [(0, (0, 1))] + sites + [(len(aa_seq), (len(aa_seq), len(aa_seq)))]

        windows: List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]] = []
        max_sites = len(sites)
        for i, (lhs, lrange) in enumerate(sites[:-1]):
            max_j = min(i + self.cleavage_params.miscleavage + 1, max_sites - 1)
            for j in range(i + 1, max_j + 1):
                rhs, rrange = sites[j]
                # Do not filter by length here: caller may trim leading M
                # (producing a shorter valid peptide). Final filters are applied
                # after translational modifications.
                windows.append((lhs, rhs, lrange, rrange))
        return windows

    def _generate_sliding_windows(self, aa_seq: AminoAcidSeqRecord
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Generate fixed-length windows for sliding window mode."""
        windows: List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]] = []
        max_len = min(int(self.cleavage_params.max_length), len(aa_seq))
        min_len = self.cleavage_params.min_length
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

    def _generate_archipel_windows(self, aa_seq:AminoAcidSeqRecord,
            List, cds_start: int
            ) -> List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]]:
        """Generate archipel windows based on variant islands and flanking.

        Rules mirrored from islands graph:
        - Peptide must include at least one variant AA (island).
        - Both ends must be reference AA (flanks).
        - Any consecutive reference segment between variant islands within the
          peptide must be ≤ flanking_size.
        """
        flanking_size = self.cleavage_params.flanking_size
        n = len(aa_seq)
        if n == 0:
            return []
        if cds_start is None:
            # Without frame anchor, cannot reliably map nt→AA; bail out.
            return []

        # Build variant mask over AA positions
        var_mask = [False] * n
        for vc in self.variants:
            loc = vc.location
            # Map transcript nt to AA indices relative to cds_start
            start_nt = loc.start
            end_nt = loc.end
            if end_nt <= cds_start:
                continue
            aa_start = max(0, (start_nt - cds_start) // 3)
            # Inclusive coverage: variants spanning multiple nt map across codons
            aa_end = min(n, max(aa_start, (end_nt - cds_start + 2) // 3))
            for i in range(aa_start, aa_end):
                if 0 <= i < n:
                    var_mask[i] = True

        # Find variant islands (contiguous True runs)
        islands: List[Tuple[int, int]] = []  # [start, end) in AA indices
        i = 0
        while i < n:
            if var_mask[i]:
                j = i + 1
                while j < n and var_mask[j]:
                    j += 1
                islands.append((i, j))
                i = j
            else:
                i += 1

        if not islands:
            return []

        # Split islands into clusters where intervening reference runs ≤ flanking_size
        clusters: List[List[Tuple[int, int]]] = []
        cur_cluster: List[Tuple[int, int]] = [islands[0]]
        for k in range(1, len(islands)):
            prev = islands[k - 1]
            cur = islands[k]
            # reference run length between islands
            reef_len = cur[0] - prev[1]
            if reef_len <= flanking_size:
                cur_cluster.append(cur)
            else:
                clusters.append(cur_cluster)
                cur_cluster = [cur]
        clusters.append(cur_cluster)

        windows: List[Tuple[int, int, Tuple[int, int], Tuple[int, int]]] = []
        for cluster in clusters:
            first_start, _ = cluster[0]
            _, last_end = cluster[-1]

            # Determine available flanks: count consecutive reference AA to left/right
            left_ref = 0
            idx = first_start - 1
            while idx >= 0 and not var_mask[idx] and left_ref < flanking_size:
                left_ref += 1
                idx -= 1

            right_ref = 0
            idx = last_end
            while idx < n and not var_mask[idx] and right_ref < flanking_size:
                right_ref += 1
                idx += 1

            # Both ends must end on reference AA; require ≥1 ref flank when available
            for l_ext in range(1, left_ref + 1):
                lhs = first_start - l_ext
                # ensure lhs is reference (by construction it is)
                for r_ext in range(1, right_ref + 1):
                    rhs = last_end + r_ext
                    # Ensure intervening reefs inside cluster obey flanking_size
                    ok = True
                    for a, b in zip(cluster, cluster[1:]):
                        reef_len = b[0] - a[1]
                        if reef_len > flanking_size:
                            ok = False
                            break
                    if not ok:
                        continue
                    # Respect min/max length constraints loosely here; final
                    # peptide validity checked later in caller.
                    windows.append((lhs, rhs, (lhs, lhs + 1), (rhs, rhs + 1)))

        return windows

    def get_effective_variants(self, lhs:int, rhs:int,
            lrange:Tuple[int,int], rrange:Tuple[int,int], cds_start:int,
            variants:List[VariantRecordWithCoordinate],
            variants_stop_lost:List[Tuple[bool,bool,bool]],
            variants_stop_gain:List[Tuple[bool,bool,bool]],
            variants_silent_mutation:List[Tuple[bool,bool,bool]]
            ) -> List[VariantRecordWithCoordinate]:
        """ Identify variants that meaningfully affect the current peptide window.

        This method filters the variant list to find those that have a functional
        impact on the given peptide slice. The `cds_start` parameter is critical:
        it defines the frame origin for all position-based checks.

        When multiple in-frame starts exist, `cds_start` should be the most recent
        upstream start (found via get_last_in_frame_orf), not the current loop's
        iteration point. This ensures:
        - Variants between earlier and current starts are correctly flagged
        - Stop-loss/gain indices are computed relative to the right frame origin
        - Upstream frameshift tracking works correctly

        Example impact of wrong cds_start:
        - Starts at 0, 300, 600 (frame 0); variant at 150; current loop at 600
        - Using cds_start=600: check "600 < 150" → False, variant missed
        - Using cds_start=0: check "0 < 150" → True, variant correctly identified

        Returns:
            List of variants that overlap, affect translation, or shift the reading
            frame within the peptide window [lhs, rhs).
        """
        effective_variants:List[VariantRecordWithCoordinate] = []
        offset = 0
        query = FeatureLocation(start=lhs, end=rhs)
        start_loc = FeatureLocation(start=cds_start, end=cds_start + 3)
        rf_index = cds_start % 3
        frames_shifted = 0
        upstream_indels:List[VariantRecordWithCoordinate] = []
        is_coding = self.is_coding \
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
