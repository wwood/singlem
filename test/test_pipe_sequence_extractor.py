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
import re

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.pipe_sequence_extractor import _align_proteins_to_hmm
from singlem.sequence_classes import SeqReader

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    headers_with_extras = headers + str.split('read_names nucleotides_aligned taxonomy_by_known?')
    maxDiff = None
    two_packages = '%s %s' % (
        os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
        os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    def test__align_proteins_to_hmm(self):
        with open(path_to_data +
                  '/4.12.22seqs.spkg/4.12.22seqs/singlem_package_creatorq4droc.fasta') as f:
            proteins = list(SeqReader().readfq(f))
        hmm = path_to_data+'/4.12.22seqs.spkg/4.12.22seqs/graftmgyqgXl_search.hmm'

        alignment = _align_proteins_to_hmm(proteins, hmm)
        self.assertEqual(22, len(alignment))
        a = alignment[0]
        self.assertEqual('2512564006', a.name)
        self.assertEqual(
            '-------MAKKVAGTMKLQVAAGKANPSPPVGPALGQRGINIMEFCKAFNAKTaDLEP-----GAPCPTVITYYQDKSFSMEIKTPPASYFLKKAAKV-----K--------SGSKTPSRDTVG---------TVTTKQVREIAEAKMKDLNANDIEGAMKIILGSARSMGIEVK---------',
            a.seq)
        a2 = alignment[21]
        self.assertEqual('2519103189', a2.name)
        self.assertEqual(
            '-------VAKKVDSVVKLQIPAGKANPAPPVGPALGQAGINIMGFCKEFNAQT-QDQA-----GMIIPVEITVYEDRSFTFITKTPPAAVLLKKAAGI-----E--------TASGEPNRNKVA---------TLNRDKVKEIAELKMPDLNAADVEAAMRMVEGTARSMGIVIED--------',
            a2.seq)

if __name__ == "__main__":
    unittest.main()
