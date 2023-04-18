#!/usr/bin/env python

#=======================================================================
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
#=======================================================================

import unittest
import sys
import extern
import tempfile
import os.path

from bird_tool_utils import in_tempdir

path_to_script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'supplement_metapackage')
path_to_data = os.path.join(os.path.dirname(__file__), 'data')

singlem_base_directory = os.path.join(os.path.dirname(__file__), '..', '..', '..')
singlem_bin_directory = os.path.join(singlem_base_directory, 'bin')

class Tests(unittest.TestCase):
    def test_defined_taxonomy(self):
        with in_tempdir():
            extern.run(f"PATH=$PATH:{singlem_bin_directory} PYTHONPATH=$PYTHONPATH:{singlem_base_directory} {path_to_script} metapackage --new-genome-fasta-files {path_to_data}/GCA_011373445.1_genomic.fna --input-metapackage {path_to_data}/4.11.22seqs.gpkg.spkg.smpkg/ --output-metapackage out.smpkg --new-taxonomies {path_to_data}/GCA_011373445.1_genomic.fna.taxonomy")


if __name__ == "__main__":
    unittest.main()
