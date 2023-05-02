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

from bird_tool_utils import in_tempdir

path_to_script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'supplement_metapackage')
path_to_data = os.path.join(os.path.dirname(__file__), 'data')

singlem_base_directory = os.path.join(os.path.dirname(__file__), '..', '..', '..')
singlem_bin_directory = os.path.join(singlem_base_directory, 'bin')

# TODO: Once GTDBtk can be included in conda env (as of diamond 2.1.7 likely), remove the added PATH entry
run = f"PATH=$PATH:{singlem_bin_directory}:~/e/gtdbtk-v2.2.6/bin PYTHONPATH=$PYTHONPATH:{singlem_base_directory} {path_to_script} metapackage"
singlem = f"{singlem_bin_directory}/singlem"


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
            import IPython; IPython.embed()
            cmd = f"{run} --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --new-taxonomies {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna.taxonomy"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms	CTTAAAAAGAAACTAAAAGGT------GCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.16	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__NEW_SPECIES'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

    def test_auto_taxonomy(self):
        with in_tempdir():
            cmd = f"{run} --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms	CTTAAAAAGAAACTAAAAGGT------GCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.16	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

    def test_auto_taxonomy_with_not_all_new(self):
        with in_tempdir():
            cmd = f"{run} --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna {path_to_data}/GCA_011373445.1_genomic.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms	CTTAAAAAGAAACTAAAAGGT------GCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.16	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms',
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.06	Root; d__Archaea; p__Crenarchaeota; c__Thermoprotei; o__Desulfurococcales; f__Desulfurococcaceae; g__Thermosphaera'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

    def test_auto_taxonomy_with_not_all_new_fast(self):
        with in_tempdir():
            cmd = f"{run} --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna {path_to_data}/GCA_011373445.1_genomic.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --gtdbtk-output-directory {path_to_data}/GCA_011373445.1_genomic_and_mutated93.gtdbtk_output --output-taxonomies out.taxonomy.tsv"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.fna {path_to_data}/GCA_011373445.1_genomic.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms	CTTAAAAAGAAACTAAAAGGT------GCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.16	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms',
                '4.11.22seqs	GCA_011373445.1_genomic	CTAAAAAAGAAACTAAAAGAT------GTTCATATGAGGGTTATAAAAAACACTCTAATG	1	1.06	Root; d__Archaea; p__Crenarchaeota; c__Thermoprotei; o__Desulfurococcales; f__Desulfurococcaceae; g__Thermosphaera'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            # Assert taxonomy file is correct
            with open('out.taxonomy.tsv') as f:
                taxonomy = f.read()
            expected = "genome	taxonomy\n" \
            "GCA_011373445.1_genomic.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__DRVV01 sp011373445\n" \
            "GCA_011373445.1_genomic.mutated93_ms.fna	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__GCA_011373445.1_genomic.mutated93_ms\n"
            self.assertEqual(expected, taxonomy)


if __name__ == "__main__":
    unittest.main()
