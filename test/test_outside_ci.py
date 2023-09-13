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
import extern
import re

import pytest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from bird_tool_utils import in_tempdir

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

# TODO: Once GTDBtk can be included in conda env (as of diamond 2.1.7 likely), remove the added PATH entry
singlem_bin_directory = os.path.join(singlem_base_directory, 'bin')
run = f"GTDBTK_DATA_PATH=/work/microbiome/db/gtdb/gtdb_release207_v2 PATH=$PATH:{singlem_bin_directory} {path_to_script} supplement"
singlem = f"{singlem_bin_directory}/singlem"

class Tests(unittest.TestCase):
    '''Tests which require a full metapackage installation to run, so can't be run through GitHub actions.'''
    maxDiff = None

    # Used in test_outside_ci.py too
    def assert_equal_taxonomic_profile(self, observed_string, expected_array_or_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])

        r = re.compile(r'  +')
        if isinstance(expected_array_or_string, str):
            expected_array = list(expected_array_or_string.split("\n"))
            expected_array = [r.sub("\t", line) for line in expected_array]
            expected_array = [line.split("\t") for line in expected_array]
        else:
            expected_array = expected_array_or_string

        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    @pytest.mark.skipif(os.environ.get("SINGLEM_METAPACKAGE_PATH") is None, reason="Appear to be running in CI")
    def test_condense_cli(self):
        '''Test the condense CLI.'''
        with in_tempdir():
            cmd1 = "{} pipe --genome-fasta-files {}/methanobacteria/transcriptomes/GB_GCA_000309865.1_protein.fna --archive-otu-table {}".format(
                path_to_script, path_to_data, 'mbac.json')
            extern.run(cmd1)
            cmd = "{} condense --input-archive-otu-table mbac.json -p /dev/stdout".format(path_to_script)
            observed = extern.run(cmd)
            expected = '\n'.join(
                ['sample  coverage        taxonomy',
                'GB_GCA_000309865.1_protein      1.12    Root; d__Archaea; p__Methanobacteriota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__Methanobacterium sp000309865\n'])
            self.assert_equal_taxonomic_profile(observed, expected)

    # This test takes a long time - like 1+ hours.
    @pytest.mark.skipif(os.environ.get("SINGLEM_METAPACKAGE_PATH") is None, reason="Appear to be running in CI")
    def test_supplement_with_extra_taxon_genome_lengths(self):
        with in_tempdir():
            cmd = f"{path_to_script} --skip-taxonomy-check --hmmsearch-evalue 1e-5 --no-quality-filter --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --new-taxonomies {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna.taxonomy --checkm2-quality-file ~/git/singlem/test/data/supplement/checkm2.output/quality_report.tsv"
            extern.run(cmd)

            cmd2 = f'{singlem} pipe --genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps.fna --metapackage out.smpkg/ --otu-table /dev/stdout'
            output = extern.run(cmd2)
            expected = [
                "\t".join(self.headers),
                '4.11.22seqs	GCA_011373445.1_genomic.mutated93_ms.manually_added_nongaps	CTTAAAAAGAAACTAAAAGGTGCCGGCGCTCACATGAGGGTTCTAAAAAACACTCTAATT	1	1.18	Root; d__Archaea; p__Thermoproteota; c__Bathyarchaeia; o__B26-1; f__UBA233; g__DRVV01; s__NEW_SPECIES'
            ]
            self.assertEqualOtuTable(list([line.split("\t") for line in expected]), output)

            import IPython; IPython.embed()


if __name__ == "__main__":
    import logging
    # logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    unittest.main()