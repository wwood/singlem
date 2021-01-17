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
import tempdir
import json

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.singlem_package import SingleMPackage

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

if __name__ == "__main__":
    unittest.main()
