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
from singlem.chainsaw import Chainsaw

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_chainsaw_package(self):
        with tempdir.in_tempdir():
            Chainsaw.chainsaw(
                input_singlem_package_path = os.path.join(
                    path_to_data, '4.11.22seqs.gpkg.spkg'),
                output_singlem_package_path = 'chainsaw.spkg',
                sequence_prefix = "4.11~")
            self.assertTrue(os.path.isdir("chainsaw.spkg"))
            self.assertListEqual(
                list(io.open("chainsaw.spkg/4.11.22seqs/singlem_package_creatorB2tZGH.fasta")), 
                list(io.open(os.path.join(path_to_data, "chainsaw.fasta"))))
            self.assertListEqual(
                list(io.open("chainsaw.spkg/4.11.22seqs/4.11.22seqs.gpkg.refpkg/2_seqinfo.csv")),
                list(io.open(os.path.join(path_to_data,"chainsaw_seqinfo.csv"))))

if __name__ == "__main__":
    unittest.main()
