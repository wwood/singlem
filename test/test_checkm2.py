#!/usr/bin/env python3

# =======================================================================
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
# =======================================================================

import unittest
import os.path
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path
from singlem.checkm2 import CheckM2


class Tests(unittest.TestCase):
    maxDiff = None

    def test_simple(self):
        checkm2 = CheckM2(os.path.join(path_to_data, 'checkm2_example_quality_report.tsv'))
        self.assertEqual(['NASQAN2010_155_F_bin.27'], checkm2.genomes_of_sufficient_quality(95, 5))

        self.assertEqual(['NASQAN2010_127_B_bin.3', 'NASQAN2010_155_B_bin.3', 'NASQAN2010_155_F_bin.27'],
                         checkm2.genomes_of_sufficient_quality(70, 10))


if __name__ == "__main__":
    unittest.main()
