#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import sys
import os
import unittest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from singlem.prefilter_pad import pad_aligned_sequence, window_alignment_positions
from singlem.sequence_classes import AlignedProteinSequence


class PadTests(unittest.TestCase):
    def test_simple_full_coverage(self):
        # window at columns 2,3,4 => C,D,E; flank 2 on each side
        self.assertEqual(
            'ABCDEFG',
            pad_aligned_sequence('ABCDEFG', [2, 3, 4], flank_length=2))

    def test_left_and_right_padding(self):
        self.assertEqual(
            'XXABCDEFGXX',
            pad_aligned_sequence('ABCDEFG', [2, 3, 4], flank_length=4))

    def test_inserts_within_window_removed(self):
        # lowercase 'x' is an insert column between the window columns and so is
        # not one of the window_columns; it must not appear in the output at all.
        self.assertEqual(
            'ABDEFGX',
            pad_aligned_sequence('ABDxEFG', [2, 4, 5], flank_length=2))

    def test_deletion_within_window_padded(self):
        # '-' at a window column becomes the pad character
        self.assertEqual(
            'ABCXEFG',
            pad_aligned_sequence('ABC-EFG', [2, 3, 4], flank_length=2))

    def test_returns_none_when_window_start_not_covered(self):
        self.assertIsNone(
            pad_aligned_sequence('AB-DEFG', [2, 3, 4], flank_length=2))

    def test_returns_none_when_window_end_not_covered(self):
        self.assertIsNone(
            pad_aligned_sequence('ABCD-FG', [2, 3, 4], flank_length=2))

    def test_insert_residues_in_flank_kept_and_uppercased(self):
        # leading lower-case inserts are flank context and are upper-cased
        self.assertEqual(
            'ABCDEXX',
            pad_aligned_sequence('abCDE', [2, 3, 4], flank_length=2))

    def test_default_flank_length_gives_80(self):
        left = 'A' * 40
        window = 'CDEFGHIKLMNPQRSTVWYA'  # 20 residues
        right = 'M' * 40
        aligned = left + window + right
        window_columns = list(range(40, 60))
        result = pad_aligned_sequence(aligned, window_columns)
        self.assertEqual(80, len(result))
        self.assertEqual('A' * 30, result[:30])
        self.assertEqual(window, result[30:50])
        self.assertEqual('M' * 30, result[50:])


class WindowPositionTests(unittest.TestCase):
    def test_no_inserts(self):
        # all-match alignment: match state index maps directly to column index
        seqs = [AlignedProteinSequence('s1', 'ABCDEFGHIK')]
        self.assertEqual([2, 3, 4], window_alignment_positions(seqs, 2, 3))

    def test_global_insert_column_skipped(self):
        # one sequence has a lower-case insert at column 3, making column 3 an
        # insert column for all sequences; the window should skip it.
        seqs = [
            AlignedProteinSequence('s1', 'ABC-EFGHIK'),
            AlignedProteinSequence('s2', 'ABCdEFGHIK'),
        ]
        # match states after skipping column 3: cols 0,1,2,4,5,6,7,8,9
        # window starting at match state 2, length 3 => columns 2,4,5
        self.assertEqual([2, 4, 5], window_alignment_positions(seqs, 2, 3))


if __name__ == "__main__":
    unittest.main()
