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

import unittest
import os.path
import sys
from io import StringIO
import tempfile
import pytest
import extern
import re

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.kingfisher_sra import KingfisherSra

# Make as expensive because we don't want it to run in CI.
@pytest.mark.expensive
class Tests(unittest.TestCase):
    maxDiff = None

    def assertEqualOtuTable(self, expected_array_or_string, observed_string, no_assign_taxonomy=False):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])

        r = re.compile(r'  +')
        if isinstance(expected_array_or_string, str):
            expected_array = list(expected_array_or_string.split("\n"))
            expected_array = [r.sub("\t", line) for line in expected_array]
            expected_array = [line.split("\t") for line in expected_array]
            if no_assign_taxonomy:
                expected_array2 = [expected_array[0]]
                for line in expected_array[1:]:
                    expected_array2.append(line+[''])
                expected_array = expected_array2
        else:
            expected_array = expected_array_or_string

        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))
    
    def test_unpaired(self):
        unpaired_input = '>DRR128717.1.1 1 length=300\n\
NAAGACGGCAAGGCGGCCGTGTGGGTCGGCGCCGCGGACGACCTCTGGCAGCTGGGCAAGCCGCGCGGCT\n\
TCGGCGGCCCGTGGAAGAACACGACCGTCCAGAAAGGCGCCCCCTCCGACCCCTACCTCATGACCGGCTA\n\
CGACCAGAAGACCCTCAAGCTGACGAACCACGGGGCGAAGCCGGGACGTGTCTCCGTCCCCCTCGACGTC\n\
GCCGGCGCCGGACGCACCCCGGCCCACCCCCGCTTTGCGGTCCGCGCGGGCAAGGCGTCGGTTCCCCCCC\n\
TCCAGCGCGCCCCGACGGCC\n\
>DRR128717.2.1 2 length=300\n\
NAACAAAGTCGCTAGTGAAATTCGGGCAAGAATGAAGTAAAATAGATTCACAACTTACGCCGCTGTTCGT\n\
TACGGCGAACTCGGTAAGACGACCGTCACGGCGGGCGTGTTTCTCCAAATAGGAGAAATGCGCCCGCCGC\n\
TTTTGTTTAACCTAAAAGGAGTATCCAATGCTTGACCCCGTTGTGCAATCGGCCTTTGTTTTGATCCTTG\n'

        expected = '>DRR128717.1\n\
NAAGACGGCAAGGCGGCCGTGTGGGTCGGCGCCGCGGACGACCTCTGGCAGCTGGGCAAGCCGCGCGGCT\
TCGGCGGCCCGTGGAAGAACACGACCGTCCAGAAAGGCGCCCCCTCCGACCCCTACCTCATGACCGGCTA\
CGACCAGAAGACCCTCAAGCTGACGAACCACGGGGCGAAGCCGGGACGTGTCTCCGTCCCCCTCGACGTC\
GCCGGCGCCGGACGCACCCCGGCCCACCCCCGCTTTGCGGTCCGCGCGGGCAAGGCGTCGGTTCCCCCCC\
TCCAGCGCGCCCCGACGGCC\n\
>DRR128717.2\n\
NAACAAAGTCGCTAGTGAAATTCGGGCAAGAATGAAGTAAAATAGATTCACAACTTACGCCGCTGTTCGT\
TACGGCGAACTCGGTAAGACGACCGTCACGGCGGGCGTGTTTCTCCAAATAGGAGAAATGCGCCCGCCGC\
TTTTGTTTAACCTAAAAGGAGTATCCAATGCTTGACCCCGTTGTGCAATCGGCCTTTGTTTTGATCCTTG\n'
        with tempfile.NamedTemporaryFile() as f:
            f.write(unpaired_input.encode())
            f.flush()

            with tempfile.TemporaryDirectory() as d:
                (fwd, rev) = KingfisherSra()._split_fasta(f.name, d)
                self.assertEqual(None, rev)
                with open(fwd) as ofwd:
                    self.assertEqual(expected, ofwd.read())

    def test_paired(self):
        unpaired_input = '>DRR128717.1.1 1 length=300\n\
NAAGACGGCAAGGCGGCCGTGTGGGTCGGCGCCGCGGACGACCTCTGGCAGCTGGGCAAGCCGCGCGGCT\n\
TCGGCGGCCCGTGGAAGAACACGACCGTCCAGAAAGGCGCCCCCTCCGACCCCTACCTCATGACCGGCTA\n\
CGACCAGAAGACCCTCAAGCTGACGAACCACGGGGCGAAGCCGGGACGTGTCTCCGTCCCCCTCGACGTC\n\
GCCGGCGCCGGACGCACCCCGGCCCACCCCCGCTTTGCGGTCCGCGCGGGCAAGGCGTCGGTTCCCCCCC\n\
TCCAGCGCGCCCCGACGGCC\n\
>DRR128717.1.2 2 length=300\n\
NAACAAAGTCGCTAGTGAAATTCGGGCAAGAATGAAGTAAAATAGATTCACAACTTACGCCGCTGTTCGT\n\
TACGGCGAACTCGGTAAGACGACCGTCACGGCGGGCGTGTTTCTCCAAATAGGAGAAATGCGCCCGCCGC\n\
TTTTGTTTAACCTAAAAGGAGTATCCAATGCTTGACCCCGTTGTGCAATCGGCCTTTGTTTTGATCCTTG\n'

        expected1 = '>DRR128717.1\n\
NAAGACGGCAAGGCGGCCGTGTGGGTCGGCGCCGCGGACGACCTCTGGCAGCTGGGCAAGCCGCGCGGCT\
TCGGCGGCCCGTGGAAGAACACGACCGTCCAGAAAGGCGCCCCCTCCGACCCCTACCTCATGACCGGCTA\
CGACCAGAAGACCCTCAAGCTGACGAACCACGGGGCGAAGCCGGGACGTGTCTCCGTCCCCCTCGACGTC\
GCCGGCGCCGGACGCACCCCGGCCCACCCCCGCTTTGCGGTCCGCGCGGGCAAGGCGTCGGTTCCCCCCC\
TCCAGCGCGCCCCGACGGCC\n'
        expected2 = '>DRR128717.1\n\
NAACAAAGTCGCTAGTGAAATTCGGGCAAGAATGAAGTAAAATAGATTCACAACTTACGCCGCTGTTCGT\
TACGGCGAACTCGGTAAGACGACCGTCACGGCGGGCGTGTTTCTCCAAATAGGAGAAATGCGCCCGCCGC\
TTTTGTTTAACCTAAAAGGAGTATCCAATGCTTGACCCCGTTGTGCAATCGGCCTTTGTTTTGATCCTTG\n'
        with tempfile.NamedTemporaryFile() as f:
            f.write(unpaired_input.encode())
            f.flush()

            with tempfile.TemporaryDirectory() as d:
                (fwd, rev) = KingfisherSra()._split_fasta(f.name, d)
                with open(fwd) as ofwd:
                    self.assertEqual(expected1, ofwd.read())
                with open(rev) as ofwd:
                    self.assertEqual(expected2, ofwd.read())

    def test_sra1(self):
        '''
        Run on SRR8653040.sra
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    20       26.15    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTTACAGC    1       1.31    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra_chunk1(self):
        '''
        Run on SRR8653040.sra, which has 424064 reads. This test only runs on the first 200,000 reads.
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    13       17.00    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTTACAGC    1       1.31    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond --read-chunk-number 1 --read-chunk-size 200000'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra_chunk2(self):
        '''
        Run on SRR8653040.sra, which has 424064 reads. This test only runs on the first 200,000 reads.
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    6       7.84   Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond --read-chunk-number 2 --read-chunk-size 200000'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra_chunk3(self):
        '''
        Run on SRR8653040.sra, which has 424064 reads. This test only runs on the first 200,000 reads.
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    1       1.31    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond --read-chunk-number 3 --read-chunk-size 200000'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

if __name__ == "__main__":
    unittest.main()
