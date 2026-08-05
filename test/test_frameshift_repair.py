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

from singlem.frameshift_repair import (
    parse_btop, walk_btop, repair_frameshifts, shift_position,
    resolve_ambiguous_windows)


class ParseBtopTests(unittest.TestCase):
    def test_matches_only(self):
        self.assertEqual([('match', 12)], list(parse_btop('12')))

    def test_mismatch(self):
        self.assertEqual(
            [('match', 3), ('mismatch', 'A', 'V'), ('match', 5)],
            list(parse_btop('3AV5')))

    def test_gaps(self):
        # '-' as query character is a gap in the query, as subject character a
        # gap in the subject.
        self.assertEqual(
            [('mismatch', '-', 'K'), ('mismatch', 'L', '-')],
            list(parse_btop('-KL-')))

    def test_frameshifts(self):
        self.assertEqual(
            [('match', 2), ('frameshift', '/'), ('match', 4)],
            list(parse_btop('2/-4')))
        self.assertEqual(
            [('match', 2), ('frameshift', '\\'), ('match', 4)],
            list(parse_btop('2\\-4')))

    def test_multi_digit_counts(self):
        self.assertEqual([('match', 123)], list(parse_btop('123')))


class WalkBtopTests(unittest.TestCase):
    def test_no_frameshifts(self):
        # 10 matches from nt 1 consumes 30 nt, ending at 30.
        self.assertEqual([], walk_btop(1, 30, '10'))

    def test_end_position_mismatch_returns_none(self):
        self.assertIsNone(walk_btop(1, 99, '10'))

    def test_deletion_frameshift_forward(self):
        # 2 matches (6 nt, positions 1-6), then a '/' at nt 7. A '/' shifts the
        # frame back by one, so the remaining 4 matches cover 12 nt from 6.
        frameshifts = walk_btop(1, 17, '2/-4')
        self.assertEqual([(6, 'deletion')], frameshifts)

    def test_insertion_frameshift_forward(self):
        frameshifts = walk_btop(1, 19, '2\\-4')
        self.assertEqual([(6, 'insertion')], frameshifts)

    def test_reverse_strand(self):
        # qstart > qend means the reverse strand, and the walk goes downwards.
        self.assertEqual([], walk_btop(30, 1, '10'))
        # On the reverse strand the missing base sits immediately after
        # `position` in forward coordinates, one nucleotide further along than
        # on the forward strand (contrast test_deletion_frameshift_forward).
        self.assertEqual([(24, 'deletion')], walk_btop(30, 14, '2/-4'))

    def test_reverse_strand_insertion(self):
        self.assertEqual([(23, 'insertion')], walk_btop(30, 12, '2\\-4'))

    def test_gap_in_query_consumes_no_nucleotides(self):
        # A query gap ('-K') consumes no query nucleotides, so 2+2 matches plus
        # the gap spans 12 nt.
        self.assertEqual([], walk_btop(1, 12, '2-K2'))

    def test_gap_in_subject_consumes_a_codon(self):
        self.assertEqual([], walk_btop(1, 15, '2L-2'))

    def test_unparseable_btop_returns_none(self):
        self.assertIsNone(walk_btop(1, 10, '3A'))


class RepairFrameshiftsTests(unittest.TestCase):
    def test_no_frameshifts_is_identity(self):
        self.assertEqual('ACGTACGT', repair_frameshifts('ACGTACGT', []))

    def test_deletion_inserts_ambiguous_base(self):
        # A base is missing at position 4, so one is inserted there.
        self.assertEqual('ACGTNACGT', repair_frameshifts('ACGTACGT', [(4, 'deletion')]))

    def test_insertion_removes_base(self):
        self.assertEqual('ACGTCGT', repair_frameshifts('ACGTACGT', [(4, 'insertion')]))

    def test_multiple_frameshifts_applied_right_to_left(self):
        # Both edits land where they should even though the first changes the
        # length of the sequence.
        self.assertEqual(
            'ACNGTACNGT',
            repair_frameshifts('ACGTACGT', [(2, 'deletion'), (6, 'deletion')]))

    def test_mixed_frameshifts(self):
        self.assertEqual(
            'ACNGTCGT',
            repair_frameshifts('ACGTACGT', [(2, 'deletion'), (4, 'insertion')]))

    def test_out_of_range_frameshift_ignored(self):
        self.assertEqual('ACGT', repair_frameshifts('ACGT', [(99, 'deletion')]))

    def test_restores_reading_frame(self):
        # The point of the exercise: a sequence whose length is not a multiple
        # of 3 because of a single deleted base becomes one that is.
        sequence = 'ATGAAACCCGGGTTT'
        with_deletion = sequence[:6] + sequence[7:]
        self.assertEqual(len(sequence) - 1, len(with_deletion))
        repaired = repair_frameshifts(with_deletion, [(6, 'deletion')])
        self.assertEqual(len(sequence), len(repaired))
        self.assertEqual(0, len(repaired) % 3)


class ShiftPositionTests(unittest.TestCase):
    def test_no_frameshifts(self):
        self.assertEqual(10, shift_position(10, []))

    def test_earlier_deletion_shifts_right(self):
        self.assertEqual(11, shift_position(10, [(5, 'deletion')]))

    def test_earlier_insertion_shifts_left(self):
        self.assertEqual(9, shift_position(10, [(5, 'insertion')]))

    def test_later_frameshift_does_not_shift(self):
        self.assertEqual(10, shift_position(10, [(50, 'deletion')]))


class FakeWindowSequence:
    '''Stands in for UnalignedAlignedNucleotideSequence, which needs more setup
    than these tests require.'''

    def __init__(self, aligned_sequence, name=None):
        self.aligned_sequence = aligned_sequence
        self.name = name


class ResolveAmbiguousWindowsTests(unittest.TestCase):
    def test_no_ambiguous_windows(self):
        sequences = [FakeWindowSequence('ACGTACGT')]
        self.assertEqual(0, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGTACGT', sequences[0].aligned_sequence)

    def test_resolved_from_identical_window(self):
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous, FakeWindowSequence('ACGTACGT')]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGTACGT', ambiguous.aligned_sequence)

    def test_repaired_names_restricts_resolution(self):
        # An 'N' that came from the raw read rather than frameshift repair must
        # not be filled in, even if a perfect donor is available.
        raw_n = FakeWindowSequence('ACGNACGT', name='raw_read')
        repaired = FakeWindowSequence('ACGNACGT', name='repaired_read')
        sequences = [raw_n, repaired, FakeWindowSequence('ACGTACGT', name='donor')]
        self.assertEqual(
            1, resolve_ambiguous_windows(sequences, repaired_names={'repaired_read'}))
        self.assertEqual('ACGNACGT', raw_n.aligned_sequence)
        self.assertEqual('ACGTACGT', repaired.aligned_sequence)

    def test_repaired_names_none_resolves_everything(self):
        # Without provenance information, every ambiguous window is eligible,
        # matching the previous behaviour.
        ambiguous = FakeWindowSequence('ACGNACGT', name='unknown')
        sequences = [ambiguous, FakeWindowSequence('ACGTACGT', name='donor')]
        self.assertEqual(1, resolve_ambiguous_windows(sequences, repaired_names=None))
        self.assertEqual('ACGTACGT', ambiguous.aligned_sequence)

    def test_most_abundant_donor_wins(self):
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous]
        # 'ACGTACGT' differs at one position from 'ACGAACGT'; the more abundant
        # one should supply the base.
        sequences += [FakeWindowSequence('ACGAACGT')]
        sequences += [FakeWindowSequence('ACGTACGT') for _ in range(3)]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGTACGT', ambiguous.aligned_sequence)

    def test_donor_beyond_max_divergence_not_used(self):
        # The only candidate donor differs at 5 unambiguous positions, more than
        # the default maximum of 4, so the window is left alone.
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous, FakeWindowSequence('TTTTTTGT')]
        self.assertEqual(0, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGNACGT', ambiguous.aligned_sequence)

    def test_divergence_within_limit_is_used(self):
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous, FakeWindowSequence('ACGTTCGT')]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        # Only the ambiguous position is filled in; the mismatching base of the
        # read is kept, since it is what the read actually says.
        self.assertEqual('ACGTACGT', ambiguous.aligned_sequence)

    def test_no_unambiguous_windows_available(self):
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous, FakeWindowSequence('ACGNACGT')]
        self.assertEqual(0, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGNACGT', ambiguous.aligned_sequence)

    def test_different_length_windows_are_not_donors(self):
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous, FakeWindowSequence('ACGTACG')]
        self.assertEqual(0, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGNACGT', ambiguous.aligned_sequence)

    def test_multiple_ambiguous_bases_resolved(self):
        ambiguous = FakeWindowSequence('ANGTACNT')
        sequences = [ambiguous, FakeWindowSequence('ACGTACGT')]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGTACGT', ambiguous.aligned_sequence)

    def test_all_sequences_sharing_a_window_are_resolved(self):
        ambiguous = [FakeWindowSequence('ACGNACGT') for _ in range(3)]
        sequences = ambiguous + [FakeWindowSequence('ACGTACGT')]
        self.assertEqual(3, resolve_ambiguous_windows(sequences))
        for sequence in ambiguous:
            self.assertEqual('ACGTACGT', sequence.aligned_sequence)

    def test_max_divergence_zero_requires_exact_match(self):
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous, FakeWindowSequence('ACGTTCGT')]
        self.assertEqual(0, resolve_ambiguous_windows(sequences, max_divergence=0))
        self.assertEqual('ACGNACGT', ambiguous.aligned_sequence)

    def test_abundant_donor_wins_over_closer_rare_one(self):
        # The rare donor matches at every unambiguous position; the abundant one
        # differs at two, and disagrees about the ambiguous base. Abundance takes
        # precedence, so the abundant donor supplies the base even though it is
        # further away - donor choice must not be short-cut into "nearest wins".
        window = 'ACGTACGTACGTNACGTACGTACGTACGT'
        near = window.replace('N', 'A')                 # divergence 0, rare
        far = 'TT' + window[2:].replace('N', 'C')       # divergence 2, abundant
        ambiguous = FakeWindowSequence(window)
        sequences = [ambiguous, FakeWindowSequence(near)]
        sequences += [FakeWindowSequence(far) for _ in range(10)]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        # The window keeps its own bases; only the ambiguous one is filled, from
        # the abundant donor, so it becomes 'C' rather than the rare donor's 'A'.
        self.assertEqual(window.replace('N', 'C'), ambiguous.aligned_sequence)

    def test_window_of_all_ambiguous_bases(self):
        # Every block is ambiguous, so no block can be used to look donors up.
        # The donor must still be found rather than the window silently skipped.
        donor = 'ACGTACGT'
        ambiguous = FakeWindowSequence('N' * len(donor))
        sequences = [ambiguous, FakeWindowSequence(donor)]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        self.assertEqual(donor, ambiguous.aligned_sequence)

    def test_many_ambiguous_bases_still_finds_distant_donor(self):
        # Ambiguous bases spread across the window leave few blocks usable, and
        # the donor's mismatches sit in yet another block. The donor is still
        # within max_divergence, so it must be found.
        donor = 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'
        window = list(donor)
        for position in (0, 8, 16, 24):
            window[position] = 'N'
        window[31] = 'A'    # one real mismatch, in a block of its own
        ambiguous = FakeWindowSequence(''.join(window))
        sequences = [ambiguous, FakeWindowSequence(donor)]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        self.assertEqual(donor[:31] + 'A', ambiguous.aligned_sequence)

    def test_donors_of_other_lengths_are_ignored_not_fatal(self):
        # A donor set of mixed lengths must not upset the indexing.
        ambiguous = FakeWindowSequence('ACGNACGT')
        sequences = [ambiguous,
                     FakeWindowSequence('ACGTACG'),
                     FakeWindowSequence('ACGTACGTACGT'),
                     FakeWindowSequence('ACGTACGT')]
        self.assertEqual(1, resolve_ambiguous_windows(sequences))
        self.assertEqual('ACGTACGT', ambiguous.aligned_sequence)


class ResolveAmbiguousWindowsEquivalenceTests(unittest.TestCase):
    '''The indexed donor search must return exactly what an exhaustive search
    would, since it only exists to be faster.'''

    @staticmethod
    def _exhaustive(window, abundances, max_divergence):
        best_donor = None
        best_key = None
        for (donor, abundance) in abundances.items():
            if len(donor) != len(window):
                continue
            divergence = sum(1 for (w, d) in zip(window, donor)
                             if w != 'N' and w != d)
            if divergence > max_divergence:
                continue
            key = (-abundance, divergence, donor)
            if best_key is None or key < best_key:
                best_key = key
                best_donor = donor
        if best_donor is None:
            return window
        return ''.join(d if w == 'N' else w for (w, d) in zip(window, best_donor))

    def test_matches_exhaustive_search_on_random_data(self):
        import random
        from collections import Counter

        random.seed(20260802)
        for max_divergence in (0, 1, 2, 4):
            for _ in range(20):
                length = random.choice([12, 20, 60])
                donors = [''.join(random.choice('ACGT') for _ in range(length))
                          for _ in range(30)]
                # Windows derived from the donors, so that hits are common, with
                # a varying number of ambiguous bases and extra substitutions.
                windows = []
                for _ in range(15):
                    window = list(random.choice(donors))
                    for _ in range(random.randint(1, 4)):
                        window[random.randrange(length)] = 'N'
                    for _ in range(random.randint(0, 3)):
                        position = random.randrange(length)
                        if window[position] != 'N':
                            window[position] = random.choice('ACGT')
                    windows.append(''.join(window))

                sequences = ([FakeWindowSequence(d) for d in donors] +
                             [FakeWindowSequence(w) for w in windows])
                abundances = Counter(donors)
                expected = [self._exhaustive(w, abundances, max_divergence)
                            for w in windows]

                resolve_ambiguous_windows(sequences, max_divergence=max_divergence)
                actual = [s.aligned_sequence
                          for s in sequences[len(donors):]]
                self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
