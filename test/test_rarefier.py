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

# Set the random seed first up for reproducibility for predictable randomness
import random
random.seed(110)

import unittest
import os.path
from StringIO import StringIO
import sys

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.rarefier import Rarefier
from singlem.otu_table_collection import OtuTableCollection
from singlem.otu_table import OtuTable

class Tests(unittest.TestCase):
    def test_hello_world(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','1','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','1','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        rares = Rarefier().rarefy(table_collection, 2, random_generator=PredictableRandomGenerator())
        self.assertIsInstance(rares, OtuTable)
        rares = list(rares)
        self.assertEqual(2, len(rares))
        self.assertEqual('TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT', rares[0].sequence)
        self.assertEqual('TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA', rares[1].sequence)
        self.assertEqual([1,1], [e.count for e in rares])

    def test_not_enough_samples(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','1','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','2','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        rares = Rarefier().rarefy(table_collection, 2, random_generator=PredictableRandomGenerator())
        self.assertIsInstance(rares, OtuTable)
        rares = list(rares)
        self.assertEqual(1, len(rares))
        self.assertEqual('TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT', rares[0].sequence)
        self.assertEqual([2], [e.count for e in rares])

    def test_three_to_two(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
             ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','1','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
             ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','1','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['4.11.ribosomal_protein_L10','minimal','ATACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','1','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        rares = Rarefier().rarefy(table_collection, 2, random_generator=PredictableRandomGenerator())
        self.assertIsInstance(rares, OtuTable)
        rares = list(rares)
        self.assertEqual(2, len(rares))
        self.assertEqual('TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT', rares[0].sequence)
        self.assertEqual('ATACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT', rares[1].sequence)
        self.assertEqual([1,1], [e.count for e in rares])

    def test_multiple_genes(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
             ['4.11.ribosomal_protein_L11','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
             ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','1','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['4.11.ribosomal_protein_L10','minimal','ATACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','1','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        rares = Rarefier().rarefy(table_collection, 2, random_generator=PredictableRandomGenerator())
        self.assertIsInstance(rares, OtuTable)
        rares = list(rares)
        self.assertEqual(3, len(rares))
        self.assertEqual([1,1,2], [e.count for e in rares])

    def test_nothing_returned(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy']]
        exp = "\n".join(["\t".join(x) for x in e]+[''])
        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))
        self.assertEqual(0, len(list(table_collection)))
        rares = Rarefier().rarefy(table_collection, 0)
        self.assertEqual(0, len(list(rares)))
        
    def test_using_real_generator(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
             ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        rares = Rarefier().rarefy(table_collection, 2)
        self.assertIsInstance(rares, OtuTable)
        rares = list(rares)
        self.assertEqual(1, len(rares))
        self.assertEqual(2, rares[0].count)
        
    
class PredictableRandomGenerator:
    '''Generate numbers predictably, not relying on random.random
    as that can change between python versions'''
    
    def __init__(self, start_state=4):
        self.state = start_state
        self.multiplier = 1664525
        self.increment = 1013904223
 	self.modulus = 2**32

    def next(self):
        n = (self.multiplier * self.state + self.increment) % self.modulus
        self.state = n
        return n
        
    def sample(self, choose_mees, number_to_choose):
        if len(choose_mees) < number_to_choose: raise
        index_choices = list(range(number_to_choose))
        indices = []
        while len(indices) < number_to_choose:
            chosen = index_choices[self.next() % len(index_choices)]
            indices.append(chosen)
            index_choices.pop(chosen)
        return [choose_mees[i] for i in indices]
            
        

if __name__ == "__main__":
    unittest.main()
