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


import sys, os, unittest, logging
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.otu_table_collection import *
from singlem.otu_table import *
from singlem.otu_table_entry import *

class Tests(unittest.TestCase):
    def test_exclude_distinct_duplicates(self):
        collection = OtuTableCollection()
        table = OtuTable()
        e1 = OtuTableEntry()
        e1.marker = 'gene1'
        e1.sequence = 'AAT'
        e1.sample_name = 'sample1'
        e2 = OtuTableEntry()
        e2.marker = 'gene2'
        e2.sequence = 'GGG'
        e2.sample_name = 'sample1'
        table.add([e1,e2])
        collection.otu_table_objects = [table]
        out = list(collection.excluded_duplicate_distinct_genes())
        self.assertEqual(2, len(out))

        e2.marker = 'gene1'    
        table.add([e1,e2])
        collection.otu_table_objects = [table]
        out = list(collection.excluded_duplicate_distinct_genes())
        self.assertEqual(1, len(out))
        
        
                            
if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
