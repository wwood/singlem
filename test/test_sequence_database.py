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
import logging
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
         ['ribosomal_protein_S2_rpsB_gpkg','minimal_duplicate','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','16','112.2',
          'Root; k__Bacteria; p__Firmicutes; c__Bacilli; duplicate'],
         ['ribosomal_protein_S2_rpsB_gpkg','minimal_duplicate','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','160','112.2',
          'Root; k__Bacteria; p__Firmicutes; c__Bacilli; duplicate2'],
         ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','18.8',
          'Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        otu_table = "\n".join(["\t".join(x) for x in otu_table])

        with tempdir.TempDir() as tmp:
            db_path = os.path.join(tmp, 'my.sdb')

            collection = OtuTableCollection()
            collection.add_otu_table(StringIO.StringIO(otu_table))
            SequenceDatabase.create_from_otu_table(db_path, collection)

            db2 = SequenceDatabase.acquire(db_path)
            ss1 = db2.extract_sequences_by_blast_ids(['lcl|1'])
            self.assertEqual(3, len(ss1))
            self.assertEqual('ribosomal_protein_S2_rpsB_gpkg', ss1[0].marker)
            self.assertEqual('ribosomal_protein_S2_rpsB_gpkg', ss1[1].marker)
            self.assertEqual('ribosomal_protein_S2_rpsB_gpkg', ss1[2].marker)
            self.assertEqual('Root; k__Bacteria; p__Firmicutes; c__Bacilli', ss1[0].taxonomy)
            self.assertEqual('Root; k__Bacteria; p__Firmicutes; c__Bacilli; duplicate2', ss1[1].taxonomy)
            self.assertEqual('Root; k__Bacteria; p__Firmicutes; c__Bacilli; duplicate', ss1[2].taxonomy)

            ss3 = db2.extract_sequences_by_blast_ids(['lcl|3'])
            self.assertEqual(1, len(ss3))
            s3 = ss3[0]
            self.assertEqual('ribosomal_protein_L11_rplK_gpkg', s3.marker)
            self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',s3.sequence)
            self.assertEqual('Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales', s3.taxonomy)

            self.assertEqual(os.path.join(db_path,"sequences.fasta"), db2.sequences_fasta_file)

    def test_write_dereplicated_fasta_file(self):
        input_stream = StringIO.StringIO("AAA\tabc\nAAA\tabc_\nAAT\tabc2\nAAT\tyyy\nAAT\tabd")
        output_stream = StringIO.StringIO()
        SequenceDatabase.write_dereplicated_fasta_file(input_stream, output_stream)
        self.assertEqual(
            "\n".join([
                ">lcl|1 abc~abc_",
                "AAA",
                ">lcl|2 abc2~yyy~abd",
                "AAT"])+"\n",
            output_stream.getvalue())


if __name__ == "__main__":
    unittest.main()
