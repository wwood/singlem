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
import subprocess
import os.path
import tempfile
import tempdir
from string import split
import extern
import sys
import json
from StringIO import StringIO

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.singlem_package import SingleMPackage
from singlem.taxonomy_bihash import TaxonomyBihash

class Tests(unittest.TestCase):
    maxDiff = None

    def test_hello_world(self):
        singlem_package = SingleMPackage.acquire(
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'))
        taxfile = singlem_package.graftm_package().taxtastic_taxonomy_path()

        taxonomy = """tax_id,parent_id,rank,tax_name,root,kingdom,phylum,class,order,family,genus,species
Root,Root,root,Root,Root,,,,,,,
d__Bacteria,Root,kingdom,d__Bacteria,Root,d__Bacteria,,,,,,
p__Actinobacteria,d__Bacteria,phylum,p__Actinobacteria,Root,d__Bacteria,p__Actinobacteria,,,,,
c__Actinobacteria,p__Actinobacteria,class,c__Actinobacteria,Root,d__Bacteria,p__Actinobacteria,c__Actinobacteria,,,,
o__Actinomycetales,c__Actinobacteria,order,o__Actinomycetales,Root,d__Bacteria,p__Actinobacteria,c__Actinobacteria,o__Actinomycetales,,,
o__Coriobacteriales,c__Actinobacteria,order,o__Coriobacteriales,Root,d__Bacteria,p__Actinobacteria,c__Actinobacteria,o__Coriobacteriales,,,
"""

        s = StringIO(taxonomy)
        bihash = TaxonomyBihash.parse_taxtastic_taxonomy(s)
        self.assertEqual(
            {'Root': None,
             'c__Actinobacteria': 'p__Actinobacteria',
             'd__Bacteria': 'Root',
             'o__Actinomycetales': 'c__Actinobacteria',
             'o__Coriobacteriales': 'c__Actinobacteria',
             'p__Actinobacteria': 'd__Bacteria'},
            bihash.child_to_parent)
        self.assertEqual(
            {'Root': ['d__Bacteria'],
             'c__Actinobacteria': ['o__Actinomycetales', 'o__Coriobacteriales'],
             'd__Bacteria': ['p__Actinobacteria'],
             'p__Actinobacteria': ['c__Actinobacteria']},
            bihash.parent_to_children)


if __name__ == "__main__":
    unittest.main()
