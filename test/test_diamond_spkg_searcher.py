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


if __name__ == '__main__':
    unittest.main()
