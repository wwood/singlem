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
import sys
import re
from string import split
from StringIO import StringIO

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.uc_file import UCFile

class Tests(unittest.TestCase):
    def test_hello_world(self):
        e = '''S    0    60    *    *    *    *    *    1;size=4    *
H    0    60    98.3    +    0    0    60M    0;size=2    1;size=4
C    0    6    *    *    *    *    *    1;size=4    *'''
        me = [re.split('[ \t]+', l.strip()) for l in split(e, "\n")]
        me2 = "\n".join(["\t".join([f for f in line]) for line in me])
        
        r = [r for r in UCFile(StringIO(me2))]
        self.assertEqual(2, len(r))
        self.assertEqual(['H','0;size=2','1;size=4'],
                         [r[0].record_type, r[0].query, r[0].target])
        self.assertEqual(['C','1;size=4',None],
                         [r[1].record_type, r[1].query, r[1].target])        
                            
if __name__ == "__main__":
    unittest.main()
