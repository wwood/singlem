#!/usr/bin/env python

# =======================================================================
# Author:
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
# =======================================================================

import unittest
import extern
import os.path
import pytest
import tempfile
import sys

from bird_tool_utils import in_tempdir

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.metapackage import Metapackage
from singlem.sequence_classes import SeqReader

path_to_script = 'singlem'
path_to_data = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data', 'supplement'))

singlem_base_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
singlem_bin_directory = os.path.join(singlem_base_directory, 'bin')

# In pixi environments, GTDBtk is setup correctly with data, so no need to set GTDBTK_DATA_PATH here
# gtdbtk-running tests are skipped anyway unless --run-expensive is set
run = f"singlem supplement"
singlem = "singlem"


class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    headers_with_extras = headers + str.split(
        'read_names nucleotides_aligned taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method'
    )
    maxDiff = None
    two_packages = '%s %s' % (os.path.join(path_to_data,
                                           '4.11.22seqs.gpkg.spkg'), os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    def test_defined_taxonomy(self):
        with in_tempdir():
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --new-fully-defined-taxonomies {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna.taxonomy --new-taxonomy-database-name test_supplement --new-taxonomy-database-version 1.0"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --translation-table 11 --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__NEW_SPECIES'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            # Test --new-taxonomy-database-name test_supplement --new-taxonomy-database-version 1.0"
            mpkg = Metapackage.acquire('out.smpkg')
            self.assertEqual('test_supplement', mpkg.taxonomy_database_name())
            self.assertEqual('1.0', mpkg.taxonomy_database_version())

    weebill_genome = 'binchicken_co4822.157=binchicken_co4822.157~k141_197776_18_1_1_12.fna'
    weebill_base_db = os.path.join(path_to_data, 'weebill_base.syl2db')

    def _make_weebill_input_metapackage(self, weebill_dbs=None):
        weebill_dbs = weebill_dbs or [self.weebill_base_db]
        extern.run(
            f"{singlem} metapackage --singlem-packages {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/4.11.22seqs.gpkg.spkg "
            f"--nucleotide-sdb {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/22seqs.sdb --no-taxon-genome-lengths "
            f"--weebill-db {' '.join(weebill_dbs)} --weebill-c {' '.join(['100'] * len(weebill_dbs))} "
            f"--metapackage with_weebill.smpkg")
        return 'with_weebill.smpkg'

    def _supplement_cmd(self, input_metapackage, genomes=None, taxonomy_file=None, extra=''):
        genomes = genomes or [os.path.join(path_to_data, self.weebill_genome)]
        taxonomy_file = taxonomy_file or os.path.join(path_to_data, self.weebill_genome + '.taxonomy')
        return (f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication "
                f"--skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter "
                f"--new-genome-fasta-files {' '.join(genomes)} --input-metapackage {input_metapackage} "
                f"--output-metapackage out.smpkg --new-fully-defined-taxonomies {taxonomy_file} {extra}")

    def _weebill_num_genomes(self, syl2db):
        import yaml
        extern.run(f"weebill inspect -o inspect.yaml {syl2db}")
        with open('inspect.yaml') as f:
            return yaml.safe_load(f)[0]['num_genomes']

    def test_weebill_db_in_input_metapackage_is_extended(self):
        with in_tempdir():
            mpkg = self._make_weebill_input_metapackage()
            extern.run(self._supplement_cmd(mpkg))
            out = Metapackage.acquire('out.smpkg')
            dbs = out.weebill_databases()
            self.assertEqual(1, len(dbs))
            self.assertEqual('weebill_base.syl2db', os.path.basename(dbs[0][0]))
            self.assertEqual(100, dbs[0][1])
            self.assertEqual(1, self._weebill_num_genomes(self.weebill_base_db))
            self.assertEqual(2, self._weebill_num_genomes(dbs[0][0]))

    def test_weebill_db_end_to_end_annotation(self):
        # weebill coverage of the supplemented genome maps to the taxonomy
        # supplement gave it
        from singlem.weebill import WeebillProfiler
        with in_tempdir():
            mpkg = self._make_weebill_input_metapackage()
            extern.run(self._supplement_cmd(mpkg))
            out = Metapackage.acquire('out.smpkg')
            genome_path = os.path.join(path_to_data, self.weebill_genome)
            with open('raw.tsv', 'w') as f:
                f.write('Sample_file\tGenome_file\tTrue_cov\n')
                f.write(f's1\t{genome_path}\t3.5\n')
            WeebillProfiler().annotate('raw.tsv', out, 'annotated.tsv')
            with open(os.path.join(path_to_data, self.weebill_genome + '.taxonomy')) as f:
                expected_taxonomy = f.read().strip().split('\t')[1]
            with open('annotated.tsv') as f:
                lines = f.read().splitlines()
            self.assertEqual('Sample_file\ttaxonomy\tTrue_cov', lines[0])
            self.assertEqual(2, len(lines))
            sample, taxonomy, coverage = lines[1].split('\t')
            self.assertEqual(('s1', '3.5'), (sample, coverage))
            self.assertEqual([t.strip() for t in expected_taxonomy.split(';')],
                             [t.strip() for t in taxonomy.split(';')])

    def test_weebill_genome_without_marker_hits_not_added(self):
        with in_tempdir():
            mpkg = self._make_weebill_input_metapackage()
            no_marker = os.path.join(path_to_data, 'no_marker_genome.fna')
            with open(os.path.join(path_to_data, self.weebill_genome + '.taxonomy')) as f:
                taxonomy = f.read().strip().split('\t')[1]
            with open('taxonomy.tsv', 'w') as f:
                f.write(f'{self.weebill_genome}\t{taxonomy}\n')
                f.write(f'no_marker_genome.fna\t{taxonomy.replace("s__", "s__NO_MARKER_")}\n')
            extern.run(self._supplement_cmd(
                mpkg, genomes=[os.path.join(path_to_data, self.weebill_genome), no_marker],
                taxonomy_file='taxonomy.tsv'))
            dbs = Metapackage.acquire('out.smpkg').weebill_databases()
            self.assertEqual(2, self._weebill_num_genomes(dbs[0][0]))

    def test_no_weebill(self):
        with in_tempdir():
            mpkg = self._make_weebill_input_metapackage()
            extern.run(self._supplement_cmd(mpkg, extra='--no-weebill'))
            self.assertEqual([], Metapackage.acquire('out.smpkg').weebill_databases())

    def test_weebill_db_base(self):
        with in_tempdir():
            extern.run(self._supplement_cmd(
                f'{path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/',
                extra=f'--weebill-db-base {self.weebill_base_db} --weebill-c 100'))
            out = Metapackage.acquire('out.smpkg')
            self.assertEqual(7, out.version)
            dbs = out.weebill_databases()
            self.assertEqual(1, len(dbs))
            self.assertEqual(100, dbs[0][1])
            self.assertEqual(2, self._weebill_num_genomes(dbs[0][0]))

    def test_weebill_db_base_with_bundled_db_croaks(self):
        with in_tempdir():
            mpkg = self._make_weebill_input_metapackage()
            with self.assertRaises(extern.ExternCalledProcessError):
                extern.run(self._supplement_cmd(
                    mpkg, extra=f'--weebill-db-base {self.weebill_base_db} --weebill-c 100'))

    def test_two_bundled_weebill_dbs_croaks(self):
        with in_tempdir():
            dummy = os.path.join(os.path.dirname(__file__), 'data', 'dummy.syl2db')
            mpkg = self._make_weebill_input_metapackage([self.weebill_base_db, dummy])
            with self.assertRaises(extern.ExternCalledProcessError):
                extern.run(self._supplement_cmd(mpkg))

    def test_weebill_genome_already_in_db_croaks(self):
        # Adding a genome whose file name is already in the weebill database
        # fails rather than adding it twice
        with in_tempdir():
            input_mpkg = f'{path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/'
            extern.run(self._supplement_cmd(
                input_mpkg, extra=f'--weebill-db-base {self.weebill_base_db} --weebill-c 100'))
            os.rename('out.smpkg', 'first.smpkg')
            grown = Metapackage.acquire('first.smpkg').weebill_databases()[0][0]
            with self.assertRaises(extern.ExternCalledProcessError) as cm:
                extern.run(self._supplement_cmd(
                    input_mpkg, extra=f'--weebill-db-base {grown} --weebill-c 100'))
            self.assertIn('already added by an earlier supplement', str(cm.exception))

    def test_output_matched_protein_sequences(self):
        with in_tempdir():
            gene_definitions = os.path.join(os.getcwd(), 'gene_definitions.tsv')
            genome_fasta = os.path.join(path_to_data, 'GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna')
            transcript_fasta = os.path.join(path_to_data, 'GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna.ffn')
            protein_fasta = os.path.join(path_to_data, 'GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna.faa')
            with open(gene_definitions, 'w') as f:
                f.write("genome_fasta\ttranscript_fasta\tprotein_fasta\n")
                f.write(f"{genome_fasta}\t{transcript_fasta}\t{protein_fasta}\n")

            matched_proteins = os.path.join(os.getcwd(), 'matched_proteins.faa')
            cmd = (
                f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths "
                f"--no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter "
                f"--new-genome-fasta-files {genome_fasta} --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ "
                f"--output-metapackage out.smpkg --new-fully-defined-taxonomies "
                f"{path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna.taxonomy "
                f"--new-taxonomy-database-name test_supplement --new-taxonomy-database-version 1.0 "
                f"--gene-definitions {gene_definitions} --output-matched-protein-sequences {matched_proteins}"
            )
            extern.run(cmd)

            mpkg = Metapackage.acquire(f"{path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/")
            package_basenames = {spkg.graftm_package_basename() for spkg in mpkg.singlem_packages}

            with open(protein_fasta) as f:
                input_protein_names = {name for name, _, _ in SeqReader().readfq(f)}

            matched_records = []
            with open(matched_proteins) as f:
                matched_records = list(SeqReader().readfq(f))

            self.assertGreater(len(matched_records), 0)
            for name, _, _ in matched_records:
                genome_name, protein_name, package_basename = name.split('‡', 2)
                self.assertEqual(
                    genome_name,
                    'GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps',
                )
                self.assertIn(protein_name, input_protein_names)
                self.assertIn(package_basename, package_basenames)

    @pytest.mark.expensive
    def test_auto_taxonomy(self):
        with in_tempdir():
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

    @pytest.mark.expensive
    def test_auto_taxonomy_with_not_all_new(self):
        with in_tempdir():
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms',
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.00	Root; d__Archaea; p__Crenarchaeota; c__Thermoprotei; o__Desulfurococcales; f__Desulfurococcaceae; g__Thermosphaera'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

    def test_auto_taxonomy_with_not_all_new_fast(self):
        with in_tempdir():
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --gtdbtk-output-directory {path_to_data}/GCA_011373445.1_genomic_and_mutated93.manually_added_nongaps.gtdbtk_output --output-taxonomies out.taxonomy.tsv"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --translation-table 11 --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps',
                # The row below should really be to species level, but there's
                # differences in the alignment from transcripts and genome
                # input, where an AA moves from one side of the gap to another.
                # To fix this could shunt the gaps to the back, but not
                # implemented at the moment as low priority.
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            # Assert taxonomy file is correct
            with open('out.taxonomy.tsv') as f:
                taxonomy = f.read()
            expected = "genome	taxonomy\n" \
                "GCA_011373445.1_genomic.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__DRVV01 sp011373445\n" \
                "GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps\n"
            self.assertEqual(expected, taxonomy)

    def test_checkm2_quality_filter(self):
        with in_tempdir():
            # Modify the gtdbtk output so genome isn't assigned to species level, otherwise checkm filter removes the only novel genome
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --gtdbtk-output-directory {path_to_data}/GCA_011373445.1_genomic_and_mutated93.manually_added_nongaps.gtdbtk_output.manually_no_species --output-taxonomies out.taxonomy.tsv --checkm2-quality-file {path_to_data}/checkm2.output/quality_report.tsv --skip-taxonomy-check"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --translation-table 11 --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01',
                # The row below should really be to species level, but there's
                # differences in the alignment from transcripts and genome
                # input, where an AA moves from one side of the gap to another.
                # To fix this could shunt the gaps to the back, but not
                # implemented at the moment as low priority.
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            # Assert taxonomy file is correct
            with open('out.taxonomy.tsv') as f:
                taxonomy = f.read()
            expected = "genome	taxonomy\n" \
                "GCA_011373445.1_genomic.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic\n"
            self.assertEqual(expected, taxonomy)

            # Test lack of --new-taxonomy-database-name test_supplement --new-taxonomy-database-version 1.0"
            mpkg = Metapackage.acquire('out.smpkg')
            self.assertEqual('custom_taxonomy_database', mpkg.taxonomy_database_name())
            self.assertEqual(None, mpkg.taxonomy_database_version())

    def test_taxonomy_file_with_not_all_new_fast(self):
        with in_tempdir():
            with open('new_taxonomies', 'w') as f:
                f.write('\t'.join(['GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna',
                    'd__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__\n']))
                f.write('\t'.join(['GCA_011373445.1_genomic.fna',
                    'd__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__DRVV01 sp011373445\n']))
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --taxonomy-file new_taxonomies --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.fna {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --output-taxonomies out.taxonomy.tsv --skip-taxonomy-check --threads 2"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --translation-table 11 --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps',
                # The row below should really be to species level, but there's
                # differences in the alignment from transcripts and genome
                # input, where an AA moves from one side of the gap to another.
                # To fix this could shunt the gaps to the back, but not
                # implemented at the moment as low priority.
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            # Assert taxonomy file is correct
            with open('out.taxonomy.tsv') as f:
                taxonomy = f.read()
            expected = "genome	taxonomy\n" \
                "GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps\n" \
                "GCA_011373445.1_genomic.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__DRVV01 sp011373445\n"
            self.assertEqual(expected, taxonomy)

    def test_provided_gene_calls(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with open(os.path.join(tempdir, 'new_taxonomies'), 'w') as f:
                f.write('\t'.join(['GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna',
                    'd__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__\n']))
                f.write('\t'.join(['GCA_011373445.1_genomic.fna',
                    'd__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__DRVV01 sp011373445\n']))
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --taxonomy-file {tempdir}/new_taxonomies --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.fna {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage {tempdir}/out.smpkg --output-taxonomies {tempdir}/out.taxonomy.tsv --skip-taxonomy-check --threads 2 --gene-definitions {path_to_data}/gene_calls.tsv"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --translation-table 11 --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage {tempdir}/out.smpkg --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps',
                # The row below should really be to species level, but there's
                # differences in the alignment from transcripts and genome
                # input, where an AA moves from one side of the gap to another.
                # To fix this could shunt the gaps to the back, but not
                # implemented at the moment as low priority.
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            # Assert taxonomy file is correct
            with open(f'{tempdir}/out.taxonomy.tsv') as f:
                taxonomy = f.read()
            expected = "genome	taxonomy\n" \
                "GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps\n" \
                "GCA_011373445.1_genomic.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__DRVV01 sp011373445\n"
            self.assertEqual(expected, taxonomy)

    def test_tilde_in_genome_name(self):
        with in_tempdir():
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/binchicken_co4822.157=binchicken_co4822.157~k141_197776_18_1_1_12.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --new-fully-defined-taxonomies {path_to_data}/binchicken_co4822.157=binchicken_co4822.157~k141_197776_18_1_1_12.fna.taxonomy --new-taxonomy-database-name test_supplement --new-taxonomy-database-version 1.0"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --translation-table 11 --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.00	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__NEW_SPECIES'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

    def test_supplemented_genome_key_to_taxonomy(self):
        # weebill reports genomes by FASTA path, and these must map back to the
        # taxonomy supplement gave them even when their name is not an accession.
        from singlem.weebill import WeebillProfiler
        genome = 'binchicken_co4822.157=binchicken_co4822.157~k141_197776_18_1_1_12.fna'
        with in_tempdir():
            cmd = f"{run} --ignore-taxonomy-database-incompatibility --no-taxon-genome-lengths --no-dereplication --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/{genome} --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --new-fully-defined-taxonomies {path_to_data}/{genome}.taxonomy --new-taxonomy-database-name test_supplement --new-taxonomy-database-version 1.0"
            extern.run(cmd)
            mpkg = Metapackage.acquire('out.smpkg')
            key = WeebillProfiler()._genome_key(os.path.join('/somewhere', genome))
            with open(os.path.join(path_to_data, genome + '.taxonomy')) as f:
                expected_taxonomy = f.read().strip().split('\t')[1]
            observed = mpkg.genome_key_to_taxonomy([key])
            self.assertEqual([key], list(observed.keys()))
            self.assertEqual(
                [t.strip() for t in expected_taxonomy.split(';')],
                [t.strip() for t in observed[key].split(';')])


if __name__ == "__main__":
    unittest.main()
