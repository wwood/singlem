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

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class Tests(unittest.TestCase):
    maxDiff = None

    def test_seqs_dna(self):
        aln = '''>s1
ga-------------TATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGACGCTGATGTG
>s2 asdas
ca---------GAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGACGCTGA----
>s3
ga-------------TATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGGGCTGATGTG-
>d4
-g----------AGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGtaACTGACGCTGATG--
'''
        expected = '''TATGGAGGAACACCAGTGGC
TATGGAGGAACACCAGTGGC
TATGGAGGAACACCAGTGGC
TATGGAGGAACACCAGTGGC
'''
        with tempfile.NamedTemporaryFile() as a:
            a.write(aln)
            a.flush()
            with tempfile.NamedTemporaryFile() as n:
                n.write(aln.replace('-',''))
                n.flush()
                
                cmd = "%s --debug seqs --alignment %s --alignment_type dna --reads %s" % (
                    path_to_script, a.name, n.name)
                self.assertEqual(sorted(expected), sorted(extern.run(cmd)))

if __name__ == "__main__":
    unittest.main()
