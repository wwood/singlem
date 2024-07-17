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
import subprocess
import os.path
import tempfile
import extern
import sys
import json

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

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
        with tempfile.NamedTemporaryFile(mode='w') as a:
            a.write(aln)
            a.flush()
            with tempfile.NamedTemporaryFile() as stderr:
                cmd = "%s seqs --debug --alignment %s --alignment-type dna"\
                      " --window-size 20 2>%s" % (
                          path_to_script, a.name, stderr.name)
                stdout = extern.run(cmd)
                # This includes ignored columns at the front, which were messing things up.
                with open(stderr.name) as stde:
                    self.assertTrue(
                        'Found best section of the alignment starting from 14\n' in \
                        stde.read())
                self.assertEqual('14\n', stdout)

    def test_info_content(self):
        self.assertEqual(
            '849\n',
            extern.run(f'{path_to_script} seqs --alignment {path_to_data}/seqs/rdrp.aln --alignment-type aa --hmm {path_to_data}/seqs/graftmljiaceib_align.hmm'))

if __name__ == "__main__":
    unittest.main()
