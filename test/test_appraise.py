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

import unittest
import os.path
from string import split
import sys
from StringIO import StringIO

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.appraiser import Appraiser
from singlem.otu_table_collection import OtuTableCollection

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None

    def test_hello_world(self):
        metagenome_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual(7, a.num_found)
        self.assertEqual(4, a.num_not_found)
        self.assertEqual('minimal', a.metagenome_sample_name)
        
    def test_multiple_samples(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(2, len(app.appraisal_results))
        a = app.appraisal_results[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(7, a.num_found)
        self.assertEqual(0, a.num_not_found)
        a = app.appraisal_results[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual(0, a.num_found)
        self.assertEqual(4, a.num_not_found)

        
    
       

if __name__ == "__main__":
    unittest.main()
