import gzip
import os
import tempfile
import unittest

from unittest import mock

from singlem.utils import first_sequence_lengths, first_sra_spot_read_lengths


class Tests(unittest.TestCase):
    def test_wrapped_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'a.fa')
            with open(path, 'w') as f:
                f.write('>a\nACGT\nAC\n>b\nA\n>c\nACG\n')
            self.assertEqual([6, 1], first_sequence_lengths(path, 2))
            self.assertEqual([6, 1, 3], first_sequence_lengths(path, 5))

    def test_gzip_fastq(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'a.fq.gz')
            with gzip.open(path, 'wt') as f:
                f.write('@a\nACGT\n+\nIIII\n@b\nAC\n+\nII\n')
            self.assertEqual([4, 2], first_sequence_lengths(path, 5))
            self.assertEqual([4], first_sequence_lengths(path, 1))

    def test_fifo_is_not_read(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'fifo')
            os.mkfifo(path)
            self.assertIsNone(first_sequence_lengths(path))

    def test_sra_paired(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'SRR8653040.sra')
        self.assertEqual(
            [[251, 250], [251, 251], [251, 251], [251, 250], [250, 251]],
            first_sra_spot_read_lengths(path, 5))

    def test_sra_nanopore(self):
        # SRR41541091 is a small Oxford Nanopore run
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'SRR41541091.sra')
        self.assertEqual(
            [[1249], [1735], [1715], [1743], [1253]],
            first_sra_spot_read_lengths(path, 5))

    def test_sra_unpaired_long(self):
        vdb_dump_output = (
            ' READ_LEN: 5000\nREAD_TYPE: SRA_READ_TYPE_BIOLOGICAL\n\n'
            ' READ_LEN: 7000\nREAD_TYPE: SRA_READ_TYPE_BIOLOGICAL\n\n')
        with mock.patch('subprocess.run') as run:
            run.return_value.stdout = vdb_dump_output
            self.assertEqual([[5000], [7000]], first_sra_spot_read_lengths('x.sra', 5))

    def test_sra_technical_reads_ignored(self):
        vdb_dump_output = (
            ' READ_LEN: 8, 150\nREAD_TYPE: SRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL\n\n')
        with mock.patch('subprocess.run') as run:
            run.return_value.stdout = vdb_dump_output
            self.assertEqual([[150]], first_sra_spot_read_lengths('x.sra', 5))

    def test_sra_missing(self):
        self.assertIsNone(first_sra_spot_read_lengths('/nonexistent/x.sra', 5))


if __name__ == '__main__':
    unittest.main()
