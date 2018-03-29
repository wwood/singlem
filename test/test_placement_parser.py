
#!/usr/bin/env python

#=======================================================================
# Authors: Ben Woodcroft, Tim Lamberton.
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
import subprocess
import os.path
import tempfile
import tempdir
from string import split
from StringIO import StringIO
import extern
import sys
import json
import re

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.placement_parser import PlacementParser
from singlem.taxonomy_bihash import TaxonomyBihash

class Tests(unittest.TestCase):
    maxDiff = None

    def test_(self):
        placement = {
            "fields": [
                "classification",
                "distal_length",
                "edge_num",
                "like_weight_ratio",
                "likelihood",
                "pendant_length"
            ],
            "version": 3,
            "tree": "(((((637000311:0.36651{0},(2600255113:0.30323{1},(2513237397:0.20355{2},2561511140:0.1579{3})0.364:0.0655{4})0.617:0.05106{5})1.000:0.58532{6},(2521172697:0.61844{7},(2582581027:0.32233{8},2571042868:0.25573{9})0.949:0.15258{10})0.902:0.1263{11})0.435:0.07766{12},(2585428157:0.73027{13},(2588253501:0.29907{14},(2556921088:0.11546{15},2518645585:0.12806{16})0.960:0.14579{17})0.991:0.28541{18})0.746:0.0561{19})0.105:0.19809{20},(2576861831:0.42312{21},(2528768167:0.20363{22},2524614704:0.28321{23})0.997:0.37729{24})0.755:0.07293{25})0.545:0.12653{26},(646564583:0.4755{27},(640069326:0.63132{28},650377985:0.541{29})0.722:0.26793{30})0.983:0.79001{31},(2516653045:0.31813{32},(2523533623:0.48892{33},2579778789:0.25687{34})0.782:0.15439{35})0.934:0.60024{36}){37};",
            "placements": [
                {
                    "p": [
                        [
                            "p__Firmicutes",
                            0.289785332337,
                            24,
                            0.825736915661,
                            -1119.74125853,
                            0.102097663521
                        ],
                        [
                            "p__Firmicutes",
                            0.389679026464,
                            21,
                            0.0934323316861,
                            -1121.9202973,
                            0.103482367825
                        ],
                        [
                            "d__Bacteria",
                            8.90258789062e-06,
                            25,
                            0.0808307526527,
                            -1122.06517725,
                            0.112723096201
                        ]
                    ],
                    "nm": [
                        [
                            "HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482_2_2_1",
                            1
                        ]
                    ]
                },
                {
                    "p": [
                        [
                            "p__Firmicutes",
                            0.265430244293,
                            24,
                            0.879016471214,
                            -1107.52673141,
                            0.0889677089655
                        ],
                        [
                            "p__Firmicutes",
                            0.423113543701,
                            21,
                            0.0604937152121,
                            -1110.20299556,
                            0.10092473344
                        ],
                        [
                            "d__Bacteria",
                            8.90258789062e-06,
                            25,
                            0.0604898135741,
                            -1110.20306006,
                            0.100916068941
                        ]
                    ],
                    "nm": [
                        [
                            "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331_1_4_1",
                            1
                        ]
                    ]
                }
            ],
            "metadata": {
                "invocation": "pplacer -j 5 --verbosity 0 --out-dir /tmp/o -c test/data/4.11.22seqs.gpkg.spkg/4.11.22seqs/4.11.22seqs.gpkg.refpkg /tmp/o/combined_alignment.aln.fa"
            }
        }

        taxonomy = """tax_id,parent_id,rank,tax_name,root,kingdom,phylum,class,order,family,genus,species
Root,Root,root,Root,Root,,,,,,,
d__Bacteria,Root,kingdom,d__Bacteria,Root,d__Bacteria,,,,,,
p__Firmicutes,d__Bacteria,phylum,p__Firmicutes,Root,d__Bacteria,p__Firmicutes,,,,,
c__Bacilli,p__Firmicutes,class,c__Bacilli,Root,d__Bacteria,p__Firmicutes,c__Bacilli,,,,
c__Clostridia,p__Firmicutes,class,c__Clostridia,Root,d__Bacteria,p__Firmicutes,c__Clostridia,,,,
o__Clostridiales,c__Clostridia,order,o__Clostridiales,Root,d__Bacteria,p__Firmicutes,c__Clostridia,o__Clostridiales,,,
"""
        bihash = TaxonomyBihash.parse_taxtastic_taxonomy(StringIO(taxonomy))
        parser = PlacementParser(placement, bihash, 0.5)
        self.assertEqual(
            ["Root",'d__Bacteria','p__Firmicutes'],
            parser.otu_placement([
                'HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331_1_4_1',
                'HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482_2_2_1',
                ]))
        # Higher threshold
        parser = PlacementParser(placement, bihash, 0.95)
        self.assertEqual(
            ["Root",'d__Bacteria'],
            parser.otu_placement([
                'HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331_1_4_1',
                'HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482_2_2_1',
                ]))
        # call Firmicutes when it isn't explicitly stated in the jplace file
        placement['placements'][0]['p'][0][0] = 'o__Clostridiales'
        placement['placements'][0]['p'][1][0] = 'o__Clostridiales'
        placement['placements'][1]['p'][0][0] = 'c__Bacilli'
        placement['placements'][1]['p'][1][0] = 'c__Bacilli'
        parser = PlacementParser(placement, bihash, 0.5)
        self.assertEqual(
            ["Root",'d__Bacteria','p__Firmicutes'],
            parser.otu_placement([
                'HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331_1_4_1',
                'HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482_2_2_1',
                ]))
        self.assertEqual(
            ["Root",'d__Bacteria','p__Firmicutes','c__Bacilli'],
            parser.otu_placement([
                'HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331_1_4_1',
                ]))

if __name__ == "__main__":
    unittest.main()
