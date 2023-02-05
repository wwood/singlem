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

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class Tests(unittest.TestCase):
    output_headers = str.split('sample  bacterial_archaeal_bases        metagenome_size read_fraction')
    maxDiff = None

    def test_marine0(self):
        cmd = "{} read_fraction -p {}/read_fraction/marine0.profile  --input-metagenome-sizes {}/read_fraction/marine0.num_bases --taxon-genome-lengths-file {}/read_fraction/gtdb_mean_genome_sizes.tsv".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('\t'.join(self.output_headers)+'\n' + '\t'.join(str.split('marine0.1       16593586562.216972      17858646300.0   92.92%'))+'\n', obs)

if __name__ == "__main__":
    unittest.main()
