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


import sys, os, unittest, logging, tempfile
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.sequence_extractor import SequenceExtractor

class Tests(unittest.TestCase):
    def test_extract(self):
        fasta = '''>1
ATG
>2 comment
AAAAA
'''
        with tempfile.NamedTemporaryFile() as f:
            f.write(fasta)
            f.flush()
            with tempfile.NamedTemporaryFile() as g:
                SequenceExtractor().extract(['1'], f.name, g.name)
                self.assertEqual(['>1','ATG',''],
                                 open(g.name).read().split("\n"))
                SequenceExtractor().extract(['1','2'], f.name, g.name)
                self.assertEqual(['>1','ATG','>2 comment','AAAAA',''],
                                 open(g.name).read().split("\n"))

    def test_extract_fwd_and_revcom(self):
        fasta = '''>1
ATG
>2 comment
AAAAA
'''
        with tempfile.NamedTemporaryFile() as f:
            f.write(fasta)
            f.flush()
            with tempfile.NamedTemporaryFile() as g:
                SequenceExtractor().extract_forward_and_reverse_complement(
                    ['1'],['1','2'], f.name, g.name)
                self.assertEqual(['>1','ATG','>1','CAT','>2 comment','TTTTT',''],
                                 open(g.name).read().split("\n"))

                            
if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
