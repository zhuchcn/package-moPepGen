""" ORF tracking for brute force peptide calling """
from __future__ import annotations
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from typing import List

class ORFTracker:
    """ Tracker for all ORFs in a given transcript. """
    def __init__(self, orf_starts:List[int]):
        self.orf_starts = orf_starts or []
        self.reading_frames = [
            [], # reading frame 0
            [], # reading frame 1
            []  # reading frame 2
        ]
        self.orf_to_index = {orf:i for i,orf in enumerate(self.orf_starts)}
        self.sort_by_reading_frame()

    def filter_orfs_before(self, pos:int):
        """ Filter ORFs that start before the given position. """
        self.orf_starts = [i for i in self.orf_starts if i <= pos]
        self.reading_frames = [[], [], []]
        self.orf_to_index = {orf:i for i,orf in enumerate(self.orf_starts)}
        self.sort_by_reading_frame()

    def sort_by_reading_frame(self):
        """ Sort ORF starts by reading frame. """
        for i, orf in enumerate(self.orf_starts):
            self.reading_frames[orf % 3].append(i)

    def get_next_in_frame_orf(self, orf:int):
        """ Get the next in frame ORF. """
        i = self.orf_to_index[orf]
        reading_frame = orf % 3
        for j in self.reading_frames[reading_frame]:
            if j > i:
                return self.orf_starts[j]
        return -1

    def get_last_in_frame_orf(self, orf:int, lhs:int):
        """ Get the last in frame ORF before the given position. """
        i = self.orf_to_index[orf]
        reading_frame = orf % 3
        for j in self.reading_frames[reading_frame][::-1]:
            if j <= i:
                return -1
            if self.orf_starts[j] <= lhs:
                return self.orf_starts[j]
        return -1
