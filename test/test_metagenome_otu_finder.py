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


import sys, os, unittest
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.metagenome_otu_finder import MetagenomeOtuFinder
from singlem.sequence_classes import *

class Tests(unittest.TestCase):
    def test__nucleotide_alignment(self):
        m = MetagenomeOtuFinder()
        self.assertEqual(('AAATTT---GGG',9),\
            m._nucleotide_alignment(AlignedProteinSequence('name','AC-D'), 'AAATTTGGG', [0,1,2,3], True))
        
    def test__nucleotide_alignment_include_inserts(self):
        m = MetagenomeOtuFinder()
        self.assertEqual(('AAA---GGG',9),\
            m._nucleotide_alignment(AlignedProteinSequence('name','AC-D'), 'AAATTTGGG', [0,2,3], True))
        self.assertEqual(('AAAttt---GGG',9),\
            m._nucleotide_alignment(AlignedProteinSequence('name','AC-D'), 'AAATTTGGG', [0,2,3], True, include_inserts=True))
        self.assertEqual(('AAAtttGGG',9),\
            m._nucleotide_alignment(AlignedProteinSequence('name','AC-D'), 'AAATTTGGG', [0,3], True, include_inserts=True))

    def test__nucleotide_alignment_aligned_nucleotides(self):
        m = MetagenomeOtuFinder()
        self.assertEqual(('AAA-TG',6),\
            m._nucleotide_alignment(AlignedProteinSequence('name','AAA-TTGGG'), 'AAATTGGG', [0,1,2,3,5,6], False))

    def test_find_best_window_with_nucleotides(self):
        m = MetagenomeOtuFinder()
        seqs = [
            'gaAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            'ga-------------TATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGACGCTGATGTG',
            'ca---------GAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGACGCTGA----',
            'ga-------------TATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGGGCTGATGTG-',
            '-g----------AGATATGGA---------------------------------------------------']
        s2 = [Sequence('seq%i' % i, seq) for i, seq in enumerate(seqs)]
        unaligned = {}
        for i, seq in enumerate(seqs):
            name = 'seq%i' % i
            unaligned[name] = seq.replace('-','')
            
        obs = m.find_windowed_sequences(
            s2,
            unaligned,
            5,
            False,
            False,
            14)
        self.assertEqual(['AAAAA','ATGGA','ATGGA','ATGGA','ATGGA'],
                         [o.aligned_sequence for o in obs])

        # now without a known window
        best_position = m.find_best_window(s2, 5, False)
        obs = m.find_windowed_sequences(
            s2,
            unaligned,
            5,
            False,
            False,
            best_position)
        self.assertEqual(['AAAAA','TATGG','TATGG','TATGG','TATGG'],
                         [o.aligned_sequence for o in obs])
        
        
                            
if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
