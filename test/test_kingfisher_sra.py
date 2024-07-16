
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

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.kingfisher_sra import KingfisherSra

class Tests(unittest.TestCase):
    maxDiff = None

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
                (fwd, rev) = KingfisherSra().split_fasta(f.name, d)
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
                (fwd, rev) = KingfisherSra().split_fasta(f.name, d)
                with open(fwd) as ofwd:
                    self.assertEqual(expected1, ofwd.read())
                with open(rev) as ofwd:
                    self.assertEqual(expected2, ofwd.read())

if __name__ == "__main__":
    unittest.main()
