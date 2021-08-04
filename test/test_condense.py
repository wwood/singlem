#!/usr/bin/env python3

#=======================================================================
# Authors: Rossen Zhao
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
import tempdir
import sys
import io
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.condense import Condenser
from singlem.otu_table_collection import StreamingOtuTableCollection

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data','condense')

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_condense_small(self): #TODO add singlem packages, example input and output comparators
        with tempdir.in_tempdir():
            stream = StreamingOtuTableCollection()
            stream.add_otu_table_file(os.path.join(
                path_to_data, 'small_condense_input.csv'))
            Condenser().condense(
                input_streaming_otu_table = stream,
                singlem_packages = [
                    os.path.join(path_to_data, 'S2.1.ribosomal_protein_L2_rplB.gpkg.spkg'),
                    os.path.join(path_to_data, 'S2.10.ribosomal_protein_S7.gpkg.spkg'),
                    os.path.join(path_to_data, 'S2.11.ribosomal_protein_S10_rpsJ.gpkg.spkg')],
                trim_percent = 5,
                output_otu_table = 'small_condense_output.csv',
                krona = '/dev/null')
            with open('small_condense_output.csv') as observed:
                with open(os.path.join(path_to_data, 'small_condense_output.csv')) as expected:
                    self.assertListEqual(list(observed), list(expected))
    
    def test_condense_two_tables_cmdline(self): #TODO add singlem packages, example input and output comparators
        with tempdir.in_tempdir():
            cmd = "{} condense --input-otu-tables {}/small_condense_input_another.csv --output-otu-table out --krona a.html --singlem-packages {}/*spkg".format(
                path_to_script, path_to_data, path_to_data, path_to_data
            )
            observed = extern.run(cmd)
            
            with open('out') as observed:
                with open(os.path.join(path_to_data, 'small_condense_output_2_samples.csv')) as expected:
                    self.assertListEqual(list(observed), list(expected))

if __name__ == "__main__":
    unittest.main()