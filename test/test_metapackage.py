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

path_to_script = 'singlem'
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
                self.assertEqual('{"singlem_metapackage_version": 7, "singlem_packages": ["4.11.22seqs.v3_archaea_targetted.gpkg.spkg"], "prefilter_db_path": "prefilter.fna.dmnd", "nucleotide_sdb": null, "sqlite_db_path_key": "read_taxonomies.duckdb", "taxon_genome_lengths": null, "taxonomy_database_name": "custom_taxonomy_database", "taxonomy_database_version": null, "diamond_prefilter_performance_parameters": "--block-size 0.5 --target-indexed -c1", "diamond_taxonomy_assignment_performance_parameters": "--block-size 0.5 --target-indexed -c1", "makeidx_sensitivity_params": null, "avg_num_genes_per_species": null, "weebill_dbs": []}',
                con.read())

    def test_metapackage_create_with_sdb(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --nucleotide-sdb test/data/a.sdb --metapackage {}/a.smpkg --no-taxon-genome-lengths".format(
                path_to_script, f
            )
            extern.run(cmd)
            with open(os.path.join(f, 'a.smpkg', 'CONTENTS.json')) as con:
                self.assertEqual('{"singlem_metapackage_version": 7, "singlem_packages": ["4.11.22seqs.v3_archaea_targetted.gpkg.spkg"], "prefilter_db_path": "prefilter.fna.dmnd", "nucleotide_sdb": "a.sdb", "sqlite_db_path_key": "read_taxonomies.duckdb", "taxon_genome_lengths": null, "taxonomy_database_name": "custom_taxonomy_database", "taxonomy_database_version": null, "diamond_prefilter_performance_parameters": "--block-size 0.5 --target-indexed -c1", "diamond_taxonomy_assignment_performance_parameters": "--block-size 0.5 --target-indexed -c1", "makeidx_sensitivity_params": null, "avg_num_genes_per_species": null, "weebill_dbs": []}',
                con.read())

    def test_metapackage_create_with_weebill_db(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --no-nucleotide-sdb --no-taxon-genome-lengths --weebill-db test/data/dummy.syl2db test/data/dummy2.syl2db --weebill-c 100 100 --metapackage {}/a.smpkg".format(
                path_to_script, f
            )
            extern.run(cmd)
            with open(os.path.join(f, 'a.smpkg', 'CONTENTS.json')) as con:
                self.assertEqual('{"singlem_metapackage_version": 7, "singlem_packages": ["4.11.22seqs.v3_archaea_targetted.gpkg.spkg"], "prefilter_db_path": "prefilter.fna.dmnd", "nucleotide_sdb": null, "sqlite_db_path_key": "read_taxonomies.duckdb", "taxon_genome_lengths": null, "taxonomy_database_name": "custom_taxonomy_database", "taxonomy_database_version": null, "diamond_prefilter_performance_parameters": "--block-size 0.5 --target-indexed -c1", "diamond_taxonomy_assignment_performance_parameters": "--block-size 0.5 --target-indexed -c1", "makeidx_sensitivity_params": null, "avg_num_genes_per_species": null, "weebill_dbs": [{"db": "dummy.syl2db", "c": 100}, {"db": "dummy2.syl2db", "c": 100}]}',
                    con.read())
            mp = Metapackage.acquire(os.path.join(f, 'a.smpkg'))
            self.assertEqual(7, mp.version)
            dbs = mp.weebill_databases()
            self.assertEqual(2, len(dbs))
            self.assertTrue(dbs[0][0].endswith('a.smpkg/dummy.syl2db'))
            self.assertEqual(100, dbs[0][1])
            self.assertTrue(dbs[1][0].endswith('a.smpkg/dummy2.syl2db'))
            self.assertEqual(100, dbs[1][1])
            self.assertTrue(os.path.exists(dbs[0][0]))
            self.assertTrue(os.path.exists(dbs[1][0]))
            # A metapackage without a weebill DB exposes an empty list.
            self.assertEqual([], Metapackage(package_paths=['test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/']).weebill_databases())

    def test_metapackage_non_two_stage_weebill_db_croaks(self):
        # Only 'weebill profile --two-stage' ever reads a bundled database, so a
        # plain .syldb would be copied in and then found unusable at pipe time.
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --no-nucleotide-sdb --no-taxon-genome-lengths --weebill-db test/data/a.sdb --weebill-c 100 --metapackage {}/a.smpkg".format(
                path_to_script, f
            )
            with self.assertRaises(Exception):
                extern.run(cmd)

    def test_metapackage_weebill_dbs_differing_c_croaks(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --no-nucleotide-sdb --no-taxon-genome-lengths --weebill-db test/data/dummy.syl2db test/data/dummy2.syl2db --weebill-c 200 100 --metapackage {}/a.smpkg".format(
                path_to_script, f
            )
            with self.assertRaises(Exception):
                extern.run(cmd)

    def test_metapackage_weebill_db_without_c_croaks(self):
        with tempfile.TemporaryDirectory(prefix='singlem') as f:
            cmd = "{} metapackage --singlem-packages test/data/4.11.22seqs.v3_archaea_targetted.gpkg.spkg/ --no-nucleotide-sdb --no-taxon-genome-lengths --weebill-db test/data/dummy.syl2db --metapackage {}/a.smpkg".format(
                path_to_script, f
            )
            with self.assertRaises(Exception):
                extern.run(cmd)

    REGIME3_METAPACKAGE = '/work/microbiome/db/singlem/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb'

    def test_genome_key_from_read_name(self):
        from singlem.metapackage_read_name_store import genome_key_from_read_name, genome_key_from_genome_file
        self.assertEqual('GCA_000013625.1', genome_key_from_read_name('ribosomal_protein_L2_rplB~GB_GCA_000013625.1'))
        self.assertEqual('GCF_000744455.1', genome_key_from_read_name('S3.1~RS_GCF_000744455.1'))
        # Sequences added by supplement: '<prefix>~<genome>‡<gene>'
        self.assertEqual('AMXMAG_0092.fa', genome_key_from_read_name(
            'EIF_2_alpha~AMXMAG_0092.fa‡AMXMAG_0092_000000000003_120_1_1_9'))
        # The genome name may itself contain '~'
        self.assertEqual('binchicken_co4822.157=binchicken_co4822.157~k141_197776_18_1_1_12', genome_key_from_read_name(
            '4.11.22seqs~binchicken_co4822.157=binchicken_co4822.157~k141_197776_18_1_1_12‡k141_1_2_3'))
        # An empty spkg prefix
        self.assertEqual('g1', genome_key_from_read_name('~g1‡gene'))
        # A supplemented genome named by accession resolves to the accession on
        # both sides
        self.assertEqual('GCA_011373445.1', genome_key_from_read_name('4.11.22seqs~GCA_011373445.1_genomic.mutated‡x'))
        self.assertEqual('GCA_011373445.1', genome_key_from_genome_file('/d/GCA_011373445.1_genomic.mutated.fna'))
        self.assertIsNone(genome_key_from_read_name('2513020051'))
        self.assertIsNone(genome_key_from_read_name('marker~no_dagger'))

    def test_genome_key_to_taxonomy(self):
        if not os.path.exists(self.REGIME3_METAPACKAGE):
            self.skipTest("GTDB r232 metapackage not present")
        mp = Metapackage.acquire(self.REGIME3_METAPACKAGE)
        d = mp.genome_key_to_taxonomy(['GCF_000744455.1', 'GCF_000191585.1'])
        self.assertEqual(2, len(d))
        self.assertIn('s__Methanobacterium_B sp000744455', d['GCF_000744455.1'])
        self.assertIn('s__Methanobacterium_B lacus', d['GCF_000191585.1'])

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
