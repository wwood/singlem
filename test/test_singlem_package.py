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
import json
from bird_tool_utils import in_tempdir
import tempfile

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.singlem_package import SingleMPackage
from singlem.singlem_package import SingleMPackageVersion4

class Tests(unittest.TestCase):
    maxDiff = None

    def test_get_sequence_ids(self):
        spkg = SingleMPackage.acquire(os.path.join(path_to_data,'4.11.22seqs.length30.gpkg.spkg'))
        self.assertEqual(sorted([
            '650377985',
            '637000311',
            '2521172697',
            '2588253501',
            '646564583',
            '2518645585',
            '2556921088',
            '2513237397',
            '2523533623',
            '2600255113',
            '2516653045',
            '2571042868',
            '2528768167',
            '2524614704',
            '2585428157',
            '2576861831',
            '2582581027',
            '640069326',
            '2561511140',
            '2579778789'
        ]), sorted(spkg.get_sequence_ids()))
    
    def test_compile_version4(self):
        with in_tempdir():
            output_package_path = "testv4.spkg"
            graftm_package_path = os.path.join(path_to_data, "61_otus.v3.gpkg.length10.spkg", "61_otus.v3")
            singlem_position = 20
            window_size = 60
            target_domains = ["Bacteria"]
            gene_description = "Test gene of type version 4"

            SingleMPackageVersion4.compile(
                output_package_path,
                graftm_package_path,
                singlem_position,
                window_size,
                target_domains,
                gene_description)

            pkg = SingleMPackage.acquire(output_package_path)

            with open(pkg.contents_path()) as f:
                contents_hash = json.load(f)
            
            self.assertEqual(contents_hash[SingleMPackage.VERSION_KEY], 4)
            self.assertTrue(contents_hash[SingleMPackage.TAXONOMY_HASH_KEY].startswith("taxonomy_"))

            taxonomy_hash = pkg.taxonomy_hash()
            expected = {
            '229854' : ['k__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Legionellales', 'f__Legionellaceae', 'g__Legionella'],
            '3761685' : ['k__Bacteria', 'p__OD1'],
            '3825327' : ['k__Archaea', 'p__Crenarchaeota', 'c__MHVG'],
            '426860' : ['k__Archaea', 'p__[Parvarchaeota]', 'c__[Parvarchaea]', 'o__YLA114'],
            '4363563' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rickettsiales', 'f__mitochondria'],
            '4336814' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rickettsiales', 'f__Pelagibacteraceae'],
            '823009' : ['k__Archaea', 'p__Euryarchaeota', 'c__DSEG', 'o__ArcA07'],
            '4423155' : ['k__Bacteria', 'p__OP11', 'c__OP11-1'],
            '4455990' : ['k__Archaea', 'p__Euryarchaeota', 'c__Methanomicrobia', 'o__Methanomicrobiales', 'f__Methanomicrobiaceae', 'g__Methanoculleus'],
            '801940' : ['k__Archaea', 'p__Crenarchaeota', 'c__MHVG'],
            '4459468' : ['k__Bacteria', 'p__Proteobacteria', 'c__Deltaproteobacteria', 'o__Desulfobacterales', 'f__Desulfobulbaceae'],
            '3779572' : ['k__Bacteria', 'p__Cyanobacteria', 'c__Chloroplast', 'o__Streptophyta'],
            '4251079' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rickettsiales', 'f__mitochondria', 'g__Lardizabala'],
            '1128285' : ['k__Archaea', 'p__Euryarchaeota', 'c__Methanomicrobia', 'o__Methanosarcinales', 'f__Methanosarcinaceae', 'g__Methanosarcina'],
            '2107103' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rickettsiales', 'f__mitochondria', 'g__Pavlova'],
            '3770699' : ['k__Archaea', 'p__[Parvarchaeota]', 'c__[Parvarchaea]', 'o__WCHD3-30'],
            '4391683' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Sphingomonadales', 'f__Sphingomonadaceae'],
            '4452949' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rickettsiales', 'f__mitochondria'],
            '4363260' : ['k__Bacteria', 'p__Cyanobacteria', 'c__Chloroplast', 'o__Chlorophyta'],
            '3190878' : ['k__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rickettsiales', 'f__mitochondria'],
            '696036' : ['k__Bacteria', 'p__Tenericutes', 'c__Mollicutes', 'o__Anaeroplasmatales', 'f__Anaeroplasmataceae', 'g__Asteroleplasma'],
            '1928988' : ['k__Archaea', 'p__Crenarchaeota', 'c__MBGB']
            }
            self.assertEqual(taxonomy_hash, expected)
            

if __name__ == "__main__":
    unittest.main()
