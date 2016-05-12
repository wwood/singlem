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
import sys
import os
import tempdir
import StringIO
from singlem.otu_table_collection import OtuTableCollection

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.sequence_database import SequenceDatabase

class Tests(unittest.TestCase):

    def test_cycle(self):
        otu_table = \
        [['gene','sample','sequence','num_hits','coverage','taxonomy'],
         ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','14.4',
          'Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
         ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.2',
          'Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
         ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','18.8',
          'Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        otu_table = "\n".join(["\t".join(x) for x in otu_table])
        
        
        with tempdir.TempDir() as tmp:
            db_path = os.path.join(tmp, 'my.sdb')
            
            collection = OtuTableCollection()
            collection.add_otu_table(StringIO.StringIO(otu_table))
            SequenceDatabase.create_from_otu_table(db_path, collection)
            
            db2 = SequenceDatabase.acquire(db_path)
            s1 = db2.extract_sequence(1)
            self.assertEqual('ribosomal_protein_L11_rplK_gpkg', s1.marker)
            self.assertEqual('minimal',s1.sample_name)
            self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',s1.sequence)
            self.assertEqual(7, s1.count)
            self.assertEqual('Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales', s1.taxonomy)
            
            s3 = db2.extract_sequence(3)
            self.assertEqual('GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT',s3.sequence)
            self.assertEqual('Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus', s3.taxonomy)
            
            self.assertEqual(os.path.join(db_path,"sequences.fasta"), db2.sequences_fasta_file)
                            
if __name__ == "__main__":
    unittest.main()
