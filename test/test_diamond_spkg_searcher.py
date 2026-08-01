#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from singlem.diamond_spkg_searcher import DiamondSpkgSearcher


class FakeDiamondProcess:
    def __init__(self, stdout_lines):
        self.stdout = stdout_lines
        self.stderr = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def wait(self):
        return 0


class DiamondSpkgSearcherTests(unittest.TestCase):
    def test_malformed_diamond_output_suggests_validating_input_reads(self):
        with tempfile.TemporaryDirectory() as tempdir:
            searcher = DiamondSpkgSearcher(num_threads=1, working_directory=tempdir)
            malformed_line = 'ERR15372545.59634363 CTGCCTCAAGCTGCTACGGCGTGCTCCGCGTCGGAACTTT\n'

            with patch(
                'singlem.diamond_spkg_searcher.Popen',
                return_value=FakeDiamondProcess([malformed_line]),
            ):
                with self.assertRaisesRegex(
                    Exception,
                    'Unexpected line format.*validate the integrity and format of the input FASTQ/FASTA files',
                ):
                    searcher._prefilter(
                        diamond_database='prefilter.dmnd',
                        read_files=['reads.fastq.gz'],
                        is_reverse_reads=False,
                        performance_parameters='',
                        sample_names=['sample.fna'],
                        min_orf_length=72,
                        context_window=None,
                    )

    def _run_prefilter_on_one_hit(self, tempdir, line, repair_frameshifts):
        searcher = DiamondSpkgSearcher(num_threads=1, working_directory=tempdir)
        with patch(
            'singlem.diamond_spkg_searcher.Popen',
            return_value=FakeDiamondProcess([line]),
        ):
            results = searcher._prefilter(
                diamond_database='prefilter.dmnd',
                read_files=['reads.fastq.gz'],
                is_reverse_reads=False,
                performance_parameters='',
                sample_names=['sample.fna'],
                min_orf_length=0,
                context_window=None,
                repair_frameshifts=repair_frameshifts,
            )
        with open(results[0].query_sequences_file) as f:
            return f.read()

    def test_frameshift_repair_inserts_ambiguous_base_for_a_deletion(self):
        # 30 nt query, aligned over its whole length, with a '/' frameshift
        # meaning a base is missing from the read: 2 matches (6 nt) then the
        # frameshift, then 7 matches over the remaining 21 nt (1-5 + 6..26).
        sequence = 'ATGAAACCCGGGTTTAAACCCGGGTTTAAA'
        line = 'read1\t{}\tgene~ref1\t1\t26\t2/-7\n'.format(sequence)
        with tempfile.TemporaryDirectory() as tempdir:
            output = self._run_prefilter_on_one_hit(tempdir, line, True)
        # Only the aligned region (nt 1-26) is written, with an 'N' inserted at
        # the frameshift, so it is one base longer than the region.
        repaired = output.split('\n')[1]
        self.assertEqual('ATGAAANCCCGGGTTTAAACCCGGGTT', repaired)
        self.assertEqual(27, len(repaired))
        self.assertEqual(0, len(repaired) % 3)

    def test_frameshift_repair_removes_extra_base_for_an_insertion(self):
        sequence = 'ATGAAACCCGGGTTTAAACCCGGGTTTAAA'
        line = 'read1\t{}\tgene~ref1\t1\t28\t2\\-7\n'.format(sequence)
        with tempfile.TemporaryDirectory() as tempdir:
            output = self._run_prefilter_on_one_hit(tempdir, line, True)
        # The extra base at the frameshift is dropped from the aligned region
        # (nt 1-28), leaving a whole number of codons.
        repaired = output.split('\n')[1]
        self.assertNotIn('N', repaired)
        self.assertEqual('ATGAAACCGGGTTTAAACCCGGGTTTA', repaired)
        self.assertEqual(27, len(repaired))
        self.assertEqual(0, len(repaired) % 3)

    def test_no_frameshift_repair_when_not_requested(self):
        # Without --repair-frameshifts, DIAMOND is not asked for BTOP and the
        # sequence is written through unchanged.
        sequence = 'ATGAAACCCGGGTTTAAACCCGGGTTTAAA'
        line = 'read1\t{}\tgene~ref1\t1\t30\n'.format(sequence)
        with tempfile.TemporaryDirectory() as tempdir:
            output = self._run_prefilter_on_one_hit(tempdir, line, False)
        self.assertEqual(sequence, output.split('\n')[1])

    def test_frameshift_repair_leaves_sequence_alone_when_btop_has_no_frameshift(self):
        sequence = 'ATGAAACCCGGGTTTAAACCCGGGTTTAAA'
        line = 'read1\t{}\tgene~ref1\t1\t30\t10\n'.format(sequence)
        with tempfile.TemporaryDirectory() as tempdir:
            output = self._run_prefilter_on_one_hit(tempdir, line, True)
        self.assertEqual(sequence, output.split('\n')[1])

    def test_inconsistent_btop_does_not_change_the_sequence(self):
        # If the BTOP walk does not land on qend then the semantics assumed by
        # walk_btop do not hold, and the read must be left alone rather than
        # edited at a guessed position.
        sequence = 'ATGAAACCCGGGTTTAAACCCGGGTTTAAA'
        line = 'read1\t{}\tgene~ref1\t1\t30\t2/-99\n'.format(sequence)
        with tempfile.TemporaryDirectory() as tempdir:
            output = self._run_prefilter_on_one_hit(tempdir, line, True)
        self.assertEqual(sequence, output.split('\n')[1])


if __name__ == '__main__':
    unittest.main()
