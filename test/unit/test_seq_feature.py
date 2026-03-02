"""Unit tests for SeqFeature helpers."""
import unittest
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation


class TestMatchedLocationOffsets(unittest.TestCase):
    """Test offset behavior of MatchedLocation."""

    def test_get_ref_positions_with_offsets(self):
        """Codon and DNA coordinate helpers should include offsets."""
        loc = MatchedLocation(
            query=FeatureLocation(
                start=0, end=6, reading_frame_index=0,
                start_offset=1, end_offset=2
            ),
            ref=FeatureLocation(
                start=100, end=106, seqname='ENST1', reading_frame_index=0,
                start_offset=2, end_offset=1
            )
        )

        self.assertEqual(loc.get_ref_codon_start(), 302)
        self.assertEqual(loc.get_ref_codon_end(), 319)
        self.assertEqual(loc.get_ref_dna_start(), 303)
        self.assertEqual(loc.get_ref_dna_end(), 317)

    def test_slice_preserves_edge_query_offsets(self):
        """Slicing to both boundaries keeps query start/end offsets."""
        loc = MatchedLocation(
            query=FeatureLocation(
                start=0, end=6, reading_frame_index=0,
                start_offset=1, end_offset=2
            ),
            ref=FeatureLocation(
                start=100, end=106, seqname='ENST1', reading_frame_index=0,
                start_offset=2, end_offset=1
            )
        )

        sub = loc[0:6]
        self.assertEqual(sub.query.start, 0)
        self.assertEqual(sub.query.end, 6)
        self.assertEqual(sub.query.start_offset, 1)
        self.assertEqual(sub.query.end_offset, 2)

    def test_slice_internal_resets_query_offsets(self):
        """Internal slice should reset query offsets to zero."""
        loc = MatchedLocation(
            query=FeatureLocation(
                start=0, end=6, reading_frame_index=0,
                start_offset=1, end_offset=2
            ),
            ref=FeatureLocation(
                start=100, end=106, seqname='ENST1', reading_frame_index=0,
                start_offset=2, end_offset=1
            )
        )

        sub = loc[1:5]
        self.assertEqual(sub.query.start, 0)
        self.assertEqual(sub.query.end, 4)
        self.assertEqual(sub.query.start_offset, 0)
        self.assertEqual(sub.query.end_offset, 0)

    def test_shift_keeps_query_offsets(self):
        """Shift should move query interval and keep query offsets."""
        loc = MatchedLocation(
            query=FeatureLocation(
                start=2, end=8, reading_frame_index=1,
                start_offset=2, end_offset=1
            ),
            ref=FeatureLocation(
                start=20, end=26, seqname='ENST2', reading_frame_index=1,
                start_offset=0, end_offset=2
            )
        )

        shifted = loc.shift(3)
        self.assertEqual(shifted.query.start, 5)
        self.assertEqual(shifted.query.end, 11)
        self.assertEqual(shifted.query.start_offset, 2)
        self.assertEqual(shifted.query.end_offset, 1)
        self.assertEqual(shifted.ref.start, 20)
        self.assertEqual(shifted.ref.end, 26)

