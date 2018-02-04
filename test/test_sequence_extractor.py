#!/usr/bin/env python

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


import sys, os, unittest, tempfile
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.sequence_extractor import SequenceExtractor
from graftm.sequence_io import Sequence

class Tests(unittest.TestCase):
    def test_extract_and_read(self):
        fasta = '''>1
ATG
>2 comment
AAAAA
'''
        with tempfile.NamedTemporaryFile() as f:
            f.write(fasta)
            f.flush()
            seqs = SequenceExtractor().extract_and_read(['1'], f.name)
            self.assertEqual(1, len(seqs))
            self.assertEqual('1', seqs[0].name)
            self.assertEqual('ATG', seqs[0].seq)

            seqs = SequenceExtractor().extract_and_read(['1','2'], f.name)
            self.assertEqual(2, len(seqs))
            self.assertEqual('1', seqs[0].name)
            self.assertEqual('ATG', seqs[0].seq)
            self.assertEqual('2', seqs[1].name)
            self.assertEqual('AAAAA', seqs[1].seq)

if __name__ == "__main__":
    unittest.main()
