
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

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.orf_length_checker import OrfLengthChecker

class Tests(unittest.TestCase):

    def test_bad(self):
        with tempfile.NamedTemporaryFile() as f:
            f.write('>seq\nAAAA\n'.encode())
            f.flush()
            self.assertFalse(OrfLengthChecker.check_sequence_file_contains_an_orf(
                f.name, 72
            ))

    def test_small_good(self):
        with tempfile.NamedTemporaryFile() as f:
            f.write(('>seq\n'+'A'*100+'\n').encode())
            f.flush()
            self.assertTrue(OrfLengthChecker.check_sequence_file_contains_an_orf(
                f.name, 72
            ))

    def test_bad_then_good_much_later(self):
        with tempfile.NamedTemporaryFile() as f:
            bad = '>seq\nAAAA\n'
            good = '>seq\n'+'A'*100+'\n'
            # 1000 lines only are read
            input = bad*500+good
            f.write(input.encode())
            f.flush()
            self.assertFalse(OrfLengthChecker.check_sequence_file_contains_an_orf(
                f.name, 72
            ))

    def test_bad_then_good_just_later(self):
        with tempfile.NamedTemporaryFile() as f:
            bad = '>seq\nAAAA\n'
            good = '>seq\n'+'A'*100+'\n'
            # 1000 lines only are read
            input = bad*100+good
            f.write(input.encode())
            f.flush()
            self.assertTrue(OrfLengthChecker.check_sequence_file_contains_an_orf(
                f.name, 72
            ))

    def test_gzip_good(self):
        with tempfile.NamedTemporaryFile() as f:
            f.write(('>seq\n'+'A'*100+'\n').encode())
            f.flush()
            extern.run("gzip -k {}".format(f.name))
            self.assertTrue(OrfLengthChecker.check_sequence_file_contains_an_orf(
                f.name+'.gz', 72
            ))

if __name__ == "__main__":
    unittest.main()
