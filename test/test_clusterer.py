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
import tempfile
import extern
from StringIO import StringIO
import sys

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.clusterer import Clusterer, SampleWiseClusteredOtu
from singlem.otu_table_collection import OtuTableCollection

class Tests(unittest.TestCase):
    def test_cluster_two(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTT','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        clusters = list(Clusterer().each_cluster(table_collection, 0.5))
        self.assertEqual(1, len(clusters))
        self.assertIsInstance(clusters[0], SampleWiseClusteredOtu)
        c = clusters[0]
        self.assertEqual(6, c.count)
        self.assertEqual(9.76/4*6, c.coverage)

    def test_no_cluster(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACT','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))

        clusters = list(Clusterer().each_cluster(table_collection, 1.0))
        self.assertEqual(2, len(clusters))
        self.assertIsInstance(clusters[0], SampleWiseClusteredOtu)
        c = clusters[0]
        self.assertEqual(2, c.count)
        self.assertEqual(4.88, c.coverage)
        c = clusters[1]
        self.assertEqual(4, c.count)
        self.assertEqual(9.76, c.coverage)

    def test_cluster_across_samples_via_script(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACT','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','ma','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        with tempfile.NamedTemporaryFile(prefix='singlem_cluster') as f:
            cmd = "%s summarise --cluster --cluster_id %f --input_otu_tables %s --output_otu_table /dev/stdout" % (
                path_to_script, 58.5/60, f.name)
            for l in ["\t".join(o) for o in e]:
                f.write(l+"\n")
            f.flush()
            output = extern.run(cmd)
            out_clusters = [o.split("\t") for o in output.split("\n")]
            self.assertEqual(
                [['gene', 'sample', 'sequence', 'num_hits', 'coverage', 'taxonomy'],
                 ['4.12.ribosomal_protein_L11_rplK',
                  'ma',
                  'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACA',
                  '4',
                  '9.76',
                  'Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                 ['4.12.ribosomal_protein_L11_rplK',
                  'minimal',
                  'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACA',
                  '2',
                  '4.88',
                  'Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                 ['']],
                out_clusters)

if __name__ == "__main__":
    unittest.main()
