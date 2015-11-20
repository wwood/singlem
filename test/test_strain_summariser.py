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
from singlem.strain_summariser import StrainSummariser
from singlem.otu_table_collection import OtuTableCollection

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')
    output_headers = split('type gene sample difference_in_bp sequence num_hits coverage taxonomy')

    def test_minimal(self):
        a = [self.headers,['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ]
        table = "\n".join(["\t".join(x) for x in a]+[''])
        
        e = [self.output_headers,
             ['reference','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','0','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['strain','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','3','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['strain','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(table))
        StrainSummariser().summarise_strains(\
                        table_collection = table_collection,
                        output_table_io = output)
        self.assertEqual(exp, output.getvalue())
        
    def test_taxonomy_focus(self):
        a = [self.headers,['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Archaea; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ]
        table = "\n".join(["\t".join(x) for x in a]+[''])
        
        e = [self.output_headers,
             ['reference','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','0','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria'],
             ['strain','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','3','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.set_target_taxonomy_by_string('Root; d__Bacteria')
        table_collection.add_otu_table(StringIO(table))
        StrainSummariser().summarise_strains(\
                        table_collection = table_collection,
                        output_table_io = output)
        self.assertEqual(exp, output.getvalue())

        
    def test_multiple_genes_and_samples(self):
        a = [self.headers,['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Archaea; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             
             ['2.13.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Archaea; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.13.ribosomal_protein_L11_rplK.gpkg','minimal','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria'],
             ['2.13.ribosomal_protein_L11_rplK.gpkg','minimal','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal2','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Archaea; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal2','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal2','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Aacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ]
        table = "\n".join(["\t".join(x) for x in a]+[''])
        
        e = [self.output_headers,
             ['reference','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','0','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria'],
             ['strain','2.12.ribosomal_protein_L11_rplK.gpkg','minimal','3','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['reference','2.13.ribosomal_protein_L11_rplK.gpkg','minimal','0','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','8','17.57','Root; d__Bacteria'],
             ['strain','2.13.ribosomal_protein_L11_rplK.gpkg','minimal','3','GAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['reference','2.12.ribosomal_protein_L11_rplK.gpkg','minimal2','0','AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','18.07','Root; d__Bacteria'],
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.set_target_taxonomy_by_string('Root; d__Bacteria')
        table_collection.add_otu_table(StringIO(table))
        StrainSummariser().summarise_strains(\
                        table_collection = table_collection,
                        output_table_io = output)
        self.assertEqual(exp, output.getvalue())
                            
if __name__ == "__main__":
    unittest.main()
