import io
import os
import tempfile
import unittest
from unittest.mock import patch

from singlem.diamond_spkg_searcher import DiamondSpkgSearcher


class DummyProc:
    def __init__(self, lines):
        self.stdout = iter(lines)
        self.stderr = io.StringIO("")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return False

    def wait(self):
        return 0


class DiamondSpkgSearcherTests(unittest.TestCase):
    def _run_prefilter_with_lines(self, lines, context_window=None):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_fasta = os.path.join(tmpdir, "reads.fna")
            open(input_fasta, "w").close()

            with patch('singlem.diamond_spkg_searcher.Popen',
                       side_effect=lambda *args, **kwargs: DummyProc(lines)):
                result = DiamondSpkgSearcher(
                    1, tmpdir, context_window=context_window)._prefilter(
                        "db", [input_fasta], False, "")[0]

            with open(result.query_sequences_file) as f:
                fasta_lines = [line.strip() for line in f if line.strip()]

            return fasta_lines, result.best_hits

    def test_context_window_trims_sequences_and_names(self):
        lines = ["read1\tAAAACCCCGGGGTTTT\tgene1~abc\t5\t12\t16\n"]

        fasta_lines, best_hits = self._run_prefilter_with_lines(lines, context_window=2)

        self.assertEqual(fasta_lines[0], ">read1:3-14,16••gene1")
        self.assertEqual(fasta_lines[1], "AACCCCGGGGTT")
        self.assertEqual(best_hits["read1:3-14,16••gene1"], "gene1~abc")

    def test_default_keeps_full_sequence_and_name(self):
        lines = ["read2\tACGTACGT\tgene2~abc\t1\t8\t8\n"]

        fasta_lines, best_hits = self._run_prefilter_with_lines(lines)

        self.assertEqual(fasta_lines[0], ">read2••gene2")
        self.assertEqual(fasta_lines[1], "ACGTACGT")
        self.assertEqual(best_hits["read2••gene2"], "gene2~abc")
