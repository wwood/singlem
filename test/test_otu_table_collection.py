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


import sys, os, unittest, logging
from string import split
from StringIO import StringIO

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.otu_table_collection import *
from singlem.otu_table import *
from singlem.otu_table_entry import *

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')

    def test_exclude_distinct_duplicates(self):
        collection = OtuTableCollection()
        table = OtuTable()
        e1 = OtuTableEntry()
        e1.marker = 'gene1'
        e1.sequence = 'AAT'
        e1.sample_name = 'sample1'
        e2 = OtuTableEntry()
        e2.marker = 'gene2'
        e2.sequence = 'GGG'
        e2.sample_name = 'sample1'
        table.add([e1,e2])
        collection.otu_table_objects = [table]
        out = list(collection.excluded_duplicate_distinct_genes())
        self.assertEqual(2, len(out))

        e2.marker = 'gene1'
        table.add([e1,e2])
        collection.otu_table_objects = [table]
        out = list(collection.excluded_duplicate_distinct_genes())
        self.assertEqual(1, len(out))

    def test_collapse_coupled(self):
        metagenome_otu_table = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','minimal_1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','minimal_1','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.11.ribosomal_protein_L10','minimal_2','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','5','10.0','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae']
        ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(metagenomes))

        expecteds = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','9','19.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        expected = "\n".join(["\t".join(x) for x in expecteds])+"\n"

        out = StringIO()
        table_collection.collapse_coupled().write_to(out)
        self.assertEqual(expected, out.getvalue())



if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
