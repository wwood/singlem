#!/usr/bin/env python3

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
import sys
from io import StringIO
import tempfile
import extern

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.metapackage import Metapackage
from singlem.otu_table_collection import OtuTableCollection
from singlem.taxonomy import TaxonomyUtils

class Tests(unittest.TestCase):
    maxDiff = None

    def test_metapackage_create_on_target_fasta(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --no-nucleotide-sdb --metapackage {}/a.smpkg --no-taxon-genome-lengths".format(
                path_to_script, f
            )
            extern.run(cmd)
            with open(os.path.join(f, 'a.smpkg', 'CONTENTS.json')) as con:
                self.assertEqual('{"singlem_metapackage_version": 4, "singlem_packages": ["4.11.22seqs.v3_archaea_targetted.gpkg.spkg"], "prefilter_db_path": "prefilter.fna.dmnd", "nucleotide_sdb": null, "sqlite_db_path_key": "read_taxonomies.sqlite3", "taxon_genome_lengths": null}',
                con.read())

    def test_metapackage_create_with_sdb(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --nucleotide-sdb test/data/a.sdb --metapackage {}/a.smpkg --no-taxon-genome-lengths".format(
                path_to_script, f
            )
            extern.run(cmd)
            with open(os.path.join(f, 'a.smpkg', 'CONTENTS.json')) as con:
                self.assertEqual('{"singlem_metapackage_version": 4, "singlem_packages": ["4.11.22seqs.v3_archaea_targetted.gpkg.spkg"], "prefilter_db_path": "prefilter.fna.dmnd", "nucleotide_sdb": "a.sdb", "sqlite_db_path_key": "read_taxonomies.sqlite3", "taxon_genome_lengths": null}',
                con.read())

    def test_metapackage_read_name_store(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.gpkg.spkg --nucleotide-sdb test/data/4.11.22seqs.gpkg.spkg.smpkg/22seqs.sdb --metapackage {}/a.smpkg --no-taxon-genome-lengths".format(
                path_to_script, f
            )
            extern.run(cmd)

            mp = Metapackage.acquire(os.path.join(f, 'a.smpkg'))
            
            
            # In [14]: h[list(h.keys())[3]]
            # Out[14]:
            # ['d__Bacteria',
            # 'p__Proteobacteria',
            # 'c__Betaproteobacteria',
            # 'o__Burkholderiales',
            # 'f__Comamonadaceae',
            # 'g__Variovorax',
            # 's__Variovorax_sp._CF313']

            # In [15]: h[list(h.keys())[4]]
            # Out[15]:
            # ['d__Bacteria',
            # 'p__Firmicutes',
            # 'c__Bacilli',
            # 'o__Lactobacillales',
            # 'f__Leuconostocaceae',
            # 'g__Weissella',
            # 's__Weissella_hellenica']

            # In [16]: list(h.keys())[3]
            # Out[16]: '2513020051'

            # In [17]: list(h.keys())[4]
            # Out[17]: '2585428030'

            self.assertEqual({
                '2513020051': 
                    ['d__Bacteria',
                    'p__Proteobacteria',
                    'c__Betaproteobacteria',
                    'o__Burkholderiales',
                    'f__Comamonadaceae',
                    'g__Variovorax',
                    's__Variovorax_sp._CF313'],
                '2585428030':
                    ['d__Bacteria',
                    'p__Firmicutes',
                    'c__Bacilli',
                    'o__Lactobacillales',
                    'f__Leuconostocaceae',
                    'g__Weissella',
                    's__Weissella_hellenica']
            }, mp.get_taxonomy_of_reads(['2513020051', '2585428030']))


if __name__ == "__main__":
    unittest.main()
