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

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')

    def test_minimal(self):
        a = [self.headers,['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ]
        table = sorted(["\t".join(x) for x in a]+[''])
        
        e = [self.headers,['2.12.ribosomal_protein_L11_rplK.gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['2.11.ribosomal_protein_L10.gpkg','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        exp = sorted(["\t".join(x) for x in e]+[''])

        input = StringIO(table)
        output = StringIO()
        with tempfile.NamedTemporaryFile(prefix='singlem_summarise_test') as f:
            StrainSummariser().summarise_strains(
                    otu_table = 
        otu_table_file = kwargs.pop('otu_table')
        output_table = kwargs.pop('output_table')
            f.write(table)
            f.flush()
            cmd = "%s summarise --strain_overview_table /dev/stdout --otu_table %s" % (path_to_script,
                                                                                     path_to_data)
            self.assertEqual(exp, sorted(subprocess.check_output(cmd, shell=True).split("\n")))
       
                            
if __name__ == "__main__":
    unittest.main()
