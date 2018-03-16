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
import os.path
from string import split
import sys
from StringIO import StringIO
import tempfile
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.dereplicator import Dereplicator

class Tests(unittest.TestCase):
    def test_hello_world(self):
        dereps = Dereplicator().dereplicate(
            ('1','2','3','4'), 2,
            {
                '1': ['d1','p1'],
                '2': ['d1'],
                '3': ['d1','p2'],
                '4': ['d1','p2']
            }, ())
        for _ in range(10):
            self.assertTrue(
                ['1','2','3']==sorted(dereps) or
                ['1','2','4']==sorted(dereps))

    def test_preferred_list(self):
        dereps = Dereplicator().dereplicate(
            ('1','2','3','4'), 2,
            {
                '1': ['d1','p1'],
                '2': ['d1'],
                '3': ['d1','p2'],
                '4': ['d1','p2']
            }, ('2','4'))
        for _ in range(10):
            self.assertTrue(
                ['1','2','4']==sorted(dereps))


if __name__ == "__main__":
    unittest.main()
