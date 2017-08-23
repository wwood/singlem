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
import tempfile
import extern

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.chancer import Chancer
from singlem.otu_table_collection import OtuTableCollection
from singlem.taxonomy import Taxonomy

class Tests(unittest.TestCase):
    chance_headers = split('sample total_seqs homogeneity_index')
    headers = split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None

    def test_script_hello_world(self):
        metagenome_otu_table = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTA','5','10.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
        ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        with tempfile.NamedTemporaryFile(prefix='singlem-chance') as f:
            f.write(metagenomes)
            f.flush()

            cmd = "%s chance --otu_table %s --taxonomy 'Root'" % (
                path_to_script, f.name)
            self.assertEqual(
                "\n".join(["\t".join(x) for x in [
                    self.chance_headers,
                    ['minimal','8.0',str((4.5*5/9 + 7) / 2)]
                ]])+"\n",
                extern.run(cmd)
            )

    def test_target_taxonomy(self):
        metagenome_otu_table = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTA','5','10.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
        ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(metagenomes))
        self.assertEqual(
            ["minimal\t8.0\t4.75"],
            [str(rp) for rp in Chancer().predict_samples(
                metagenomes = table_collection,
                target_taxonomy = [])]
        )
        self.assertEqual(
            ["minimal\t9.0\t2.5"],
            [str(rp) for rp in Chancer().predict_samples(
                metagenomes = table_collection,
                target_taxonomy = Taxonomy.split_taxonomy(
                    'Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae'))]
        )

if __name__ == "__main__":
    unittest.main()
