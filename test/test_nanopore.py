#!/usr/bin/env python3

#=======================================================================
# Authors: Joshua Mitchell.
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
import extern
import sys
import re

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data', 'nanopore', 'spkgs')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

TEST_ANNOY = False
try:
    import annoy
    TEST_ANNOY = True
    print("annoy found, running relevant tests", file=sys.stderr)
except ImportError:
    print("WARNING: annoy not found, skipping relevant tests", file=sys.stderr)
    pass

class Tests(unittest.TestCase):

    path_to_fastas = 'test/data/nanopore/fastas'
    path_to_expected = 'test/data/nanopore/expected'

    base_command = """
        singlem pipe -1 {} --singlem-package {} --no-assign-taxonomy --otu-table /dev/stdout
        """

    two_packages = '%s %s' % (
        os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
        os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def assertEqualOtuTable(self, expected_array_or_string, observed_string, no_assign_taxonomy=False):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])

        r = re.compile(r'  +')
        if isinstance(expected_array_or_string, str):
            expected_array = list(expected_array_or_string.split("\n"))
            expected_array = [r.sub("\t", line) for line in expected_array]
            expected_array = [line.split("\t") for line in expected_array]
            if no_assign_taxonomy:
                expected_array2 = [expected_array[0]]
                for line in expected_array[1:]:
                    expected_array2.append(line+[''])
                expected_array = expected_array2
        else:
            expected_array = expected_array_or_string

        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    def base_test(self, test_name):
        expected = open(os.path.join(self.path_to_expected, test_name + '.tsv')).read()
        
        cmd = self.base_command.format(
            os.path.join(self.path_to_fastas, test_name + '.fasta'),
            self.two_packages,
        )

        self.assertEqualOtuTable(expected, extern.run(cmd))

    # can it find a single hit?
    def test_one_hit(self):
        test_name = 'one_hit'
        self.base_test(test_name)

    # can it deal with a single hit that is duplicated in the read?
    # should only count once
    def test_one_hit_dupped(self):
        test_name = 'one_hit_dupped'
        self.base_test(test_name)

    # can it deal with duplex reads (i.e. with a reverse complement) with a single hit?
    # should only count once
    def test_one_hit_duplex(self):
        test_name = 'one_hit_duplex'
        self.base_test(test_name)

    # can it find two hits?
    def test_two_hits(self):
        test_name = 'two_hits'
        self.base_test(test_name)

    # can it deal with two hits that are both duplicated in the read?
    # should only count each hit once
    def test_two_hits_dupped(self):
        test_name = 'two_hits_dupped'
        self.base_test(test_name)

    # can it deal with duplex reads (i.e. with a reverse complement) with multiple hits?
    # should only count each hit once
    def test_two_hits_duplex(self):
        test_name = 'two_hits_duplex'
        self.base_test(test_name)


if __name__ == "__main__":
    unittest.main()
