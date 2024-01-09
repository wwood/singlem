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

import pickle
import unittest
import os.path
import sys
from io import StringIO
import json
from bird_tool_utils import in_tempdir

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.package_creator import PackageCreator

class Tests(unittest.TestCase):
    maxDiff = None

    def test_create_protein_pkg(self):
        with in_tempdir():
            PackageCreator().create(
                input_graftm_package = os.path.join(
                    path_to_data, '4.11.22seqs.gpkg.spkg', '4.11.22seqs'),
                input_taxonomy = os.path.join(
                    path_to_data, '4.11.22seqs.gpkg.spkg', 'input_taxonomy.tsv'),
                output_singlem_package = 'protein.spkg',
                hmm_position = 76,
                window_size = 63,
                target_domains = ["Bacteria"],
                gene_description = "",
                force = False)
            self.assertTrue(os.path.isdir('protein.spkg'))
            with open('protein.spkg/CONTENTS.json') as f:
                j = json.load(f)
            self.assertEqual(75, j['singlem_hmm_position'])
            self.assertEqual(63, j['singlem_window_size'])
            self.assertTrue(j['taxonomy_hash'].startswith("taxonomy"))

            with open(os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg', "taxonomy_hash.pickle"), 'rb') as file:
                expected_hash = pickle.load(file)
            with open(os.path.join("protein.spkg", j['taxonomy_hash']), 'rb') as file:
                observed_hash = pickle.load(file)
            self.assertDictEqual(expected_hash, observed_hash)

    def test_create_nuc_pkg(self):
        with in_tempdir():
            PackageCreator().create(
                input_graftm_package = os.path.join(
                    path_to_data, '61_otus.v3.gpkg.spkg', '61_otus.v3'),
                input_taxonomy = os.path.join(
                    path_to_data, '61_otus.v3.gpkg.spkg', 'input_taxonomy.tsv'),
                output_singlem_package = 'nuc.spkg',
                hmm_position = 888,
                window_size = 57,
                target_domains = ["Bacteria"],
                gene_description = "",
                force = False)
            self.assertTrue(os.path.isdir('nuc.spkg'))
            with open('nuc.spkg/CONTENTS.json') as f:
                j = json.load(f)
            self.assertEqual(887, j['singlem_hmm_position'])
            self.assertEqual(57, j['singlem_window_size'])
            self.assertTrue(j['taxonomy_hash'].startswith("taxonomy"))

            with open(os.path.join(path_to_data, '61_otus.v3.gpkg.spkg', "taxonomy_hash.pickle"), 'rb') as file:
                expected_hash = pickle.load(file)
            with open(os.path.join("nuc.spkg", j['taxonomy_hash']), 'rb') as file:
                observed_hash = pickle.load(file)
            self.assertDictEqual(expected_hash, observed_hash)



if __name__ == "__main__":
    unittest.main()
