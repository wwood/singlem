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
import itertools

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    def test_makedb_query(self):
        otu_table = [self.headers,['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        otu_table = "\n".join(["\t".join(x) for x in otu_table])

        with tempfile.NamedTemporaryFile() as f:
            f.write(otu_table)
            f.flush()

            with tempdir.TempDir() as d:
                cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
                                                                d,
                                                                f.name)
                subprocess.check_call(cmd, shell=True)

                cmd = "%s query --query_sequence %s --db %s/db" % (
                    path_to_script,
                    'CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA', # second sequence with an extra A at the end
                                                                d)

                expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                            ['unnamed_sequence','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
                expected = ["\t".join(x) for x in expected]+['']
                self.assertEqual(expected,
                                 subprocess.check_output(cmd, shell=True).split("\n"))

    def test_query_with_otu_table(self):
        with tempfile.NamedTemporaryFile() as f:
            query = [self.headers,
                     # second sequence with an extra A at the end
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
            query = "\n".join(["\t".join(x) for x in query])
            f.write(query)
            f.flush()

            cmd = "%s query --query_otu_table %s --db %s" % (
                path_to_script,
                f.name,
                os.path.join(path_to_data,'a.sdb'))

            expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                        ['minimal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
            expected = ["\t".join(x) for x in expected]+['']
            self.assertEqual(expected,
                             subprocess.check_output(cmd, shell=True).split("\n"))

    def test_query_with_otu_table_two_samples(self):
        with tempfile.NamedTemporaryFile() as f:
            query = [self.headers,
                     # second sequence with an extra A at the end
                     ['ribosomal_protein_L11_rplK_gpkg','maximal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
                     ] # converted A to T in the middle
            query = "\n".join(["\t".join(x) for x in query])
            f.write(query)
            f.flush()

            cmd = "%s query --query_otu_table %s --db %s" % (
                path_to_script,
                f.name,
                os.path.join(path_to_data,'a.sdb'))

            expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                        ['maximal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','2','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                        ['minimal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']
                        ]
            expected = ["\t".join(x) for x in expected]+['']
            observed = subprocess.check_output(cmd, shell=True).split("\n")
            self.assertEqual(expected, observed)


    def test_query_with_otu_table_two_samples_same_sequence(self):
        with tempfile.NamedTemporaryFile() as f:
            query = [self.headers,
                     # second sequence with an extra A at the end
                     ['ribosomal_protein_L11_rplK_gpkg','maximal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli']
                     ] # converted A to T in the middle
            query = "\n".join(["\t".join(x) for x in query])
            f.write(query)
            f.flush()

            with tempdir.TempDir() as d:
                cmd = "{} makedb --db {}/sdb --otu_table {}".format(
                    path_to_script, d, f.name)
                extern.run(cmd)

                cmd = "{} query --query_otu_table {} --db {}/sdb".format(
                    path_to_script,
                    f.name,
                    d)

                expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                            ['maximal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','maximal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                            ['maximal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                            ['minimal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','maximal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                            ['minimal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                            ]
                observed = subprocess.check_output(cmd, shell=True)
                self.assertEqualOtuTable(expected, observed)

    def test_fasta_query(self):
        with tempfile.NamedTemporaryFile() as f:
            query = "\n".join([">seq1 comment",'CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC',
                               ">sseq4",       'CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC',''])
            f.write(query)
            f.flush()

            cmd = "%s query --query_fasta %s --db %s" % (
                path_to_script,
                f.name,
                os.path.join(path_to_data,'a.sdb'))

            expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                        ['seq1','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','0','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                        ['sseq4','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
            expected = ["\t".join(x) for x in expected]+['']
            observed = subprocess.check_output(cmd, shell=True).split("\n")
            self.assertEqual(expected, observed)

    def test_query_with_gaps(self):
        with tempfile.NamedTemporaryFile() as f:
            query = "\n".join([">seq1 comment",'CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC',
                               ">sseq4",       'CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC',''])
            f.write(query)
            f.flush()

            cmd = "%s query --query_fasta %s --db %s" % (
                path_to_script,
                f.name,
                os.path.join(path_to_data,'a.sdb'))

            expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                        ['seq1','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','3','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                        ['sseq4','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
            expected = ["\t".join(x) for x in expected]+['']
            observed = subprocess.check_output(cmd, shell=True).split("\n")
            self.assertEqual(expected, observed)

    def test_n_vs_gap(self):
        otu_table = [self.headers,
                     # same as query
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                     # Ns instead of gaps, plus and A instead of C at end
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAANNNTATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes']]
        otu_table = "\n".join(["\t".join(x) for x in otu_table])

        with tempfile.NamedTemporaryFile() as f:
            f.write(otu_table)
            f.flush()

            with tempdir.TempDir() as d:
                #d = '/tmp/d'
                cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
                                                                d,
                                                                f.name)
                subprocess.check_call(cmd, shell=True)

                cmd = "%s query --query_sequence %s --db %s/db" % (
                    path_to_script,
                    'CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC',
                    d)

                expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                            ['unnamed_sequence','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                            ['unnamed_sequence','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','4','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAANNNTATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes']]
                expected = ["\t".join(x) for x in expected]+['']
                self.assertEqual(expected,
                                 subprocess.check_output(cmd, shell=True).split("\n"))

    def test_duplicate_seqs_then_another(self):
        '''This tests when a two samples have the same OTU, and then there's another separate OTU after that'''
        otu_table = [self.headers,
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_L11_rplK_gpkg','maximal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                     ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        otu_table = "\n".join(["\t".join(x) for x in otu_table])

        with tempfile.NamedTemporaryFile() as f:
            f.write(otu_table)
            f.flush()

            with tempdir.TempDir() as d:
                #d = '/tmp/d'
                cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
                                                                d,
                                                                f.name)
                subprocess.check_call(cmd, shell=True)

                with tempfile.NamedTemporaryFile() as infasta:
                    infasta.write(">1_\n") # same as the first sequence
                    infasta.write("GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\n")
                    infasta.write(">2_\n") # third sequence with an replacement A at the end
                    infasta.write("CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA\n")
                    infasta.flush()

                    cmd = "%s query --query_fasta %s --db %s/db" % (
                        path_to_script,
                        infasta.name,
                        d)

                    expected = [
                        ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                        ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                        ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','9','maximal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                        ['2_','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
                    expected = ["\t".join(x) for x in expected]+['']
                    self.assertEqual(expected,
                                     subprocess.check_output(cmd, shell=True).split("\n"))


    def test_divergence0(self):
        '''Zero divergence is special because it uses the sqlite DB to query, not blast'''
        otu_table = [self.headers,
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_L11_rplK_gpkg','maximal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
        otu_table = "\n".join(["\t".join(x) for x in otu_table])



        with tempfile.NamedTemporaryFile() as f:
            f.write(otu_table)
            f.flush()

            with tempdir.TempDir() as d:

                cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
                                                                d,
                                                                f.name)
                subprocess.check_call(cmd, shell=True)

                with tempfile.NamedTemporaryFile() as infasta:
                    infasta.write(">1_\n") # same as the first sequence
                    infasta.write("GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\n")
                    infasta.write(">2_\n") # same as first
                    infasta.write("GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\n")
                    infasta.write(">3_\n") # same as third in OTU table
                    infasta.write("CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC\n")
                    infasta.flush()

                    cmd = "%s query --query_fasta %s --db %s/db --max_divergence 0" % (
                        path_to_script,
                        infasta.name,
                        d)

                    expected = [
                        ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                        ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                        ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','9','maximal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                        ['2_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                        ['2_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','9','maximal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                        ['3_','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','0','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
                    expected = ["\t".join(x) for x in expected]+['']
                    self.assertEqual(sorted(expected),
                                     sorted(subprocess.check_output(cmd, shell=True).split("\n")))

    def test_query_by_sample(self):
        expected = [
            self.headers,
            ['ribosomal_protein_L11_rplK_gpkg','m2','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.10','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['ribosomal_protein_S2_rpsB_gpkg','m2','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.40','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
            ['ribosomal_protein_S17_gpkg','m2','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGA','9','20.50','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['ribosomal_protein_L11_rplK_gpkg','m3','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.10','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['ribosomal_protein_S2_rpsB_gpkg','m3','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.40','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
            ['ribosomal_protein_S17_gpkg','m3','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','19.50','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        expected = ["\t".join(x) for x in expected]+['']

        cmd = "%s query --db %s/b.sdb --sample_names m2 m3" %(path_to_script,
                                                              path_to_data)
        self.assertEqual(expected, extern.run(cmd).split('\n'))

    def test_query_by_taxonomy(self):
        expected = [
            self.headers,
            ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.10','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','19.50','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        expected = ["\t".join(x) for x in expected]+['']

        cmd = "%s query --db %s/a.sdb --taxonomy o__Bacillales" %(path_to_script,
                                                                  path_to_data)
        self.assertEqual(expected, extern.run(cmd).split('\n'))

    def test_query_subject_otu_tables(self):
        with tempfile.NamedTemporaryFile() as f:
            query = "\n".join([">seq1 comment",'CGTCGTTGGAACCCAAAAATGAAAAAATATATCTatgTCACTGAGAGAAATGGTATTTATATC',
                               ">sseq4",       'CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC',''])
            f.write(query)
            f.flush()

            subject = [
                self.headers,
                ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.1','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.4','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','19.5','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
            with tempfile.NamedTemporaryFile() as subject_f:
                subject_f.write("\n".join(["\t".join(x) for x in subject]+['']))
                subject_f.flush()

                cmd = "%s query --query_fasta %s --subject_otu_tables %s" % (
                    path_to_script,
                    f.name,
                    subject_f.name)

                expected = [
                    ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
                    ['seq1','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTatgTCACTGAGAGAAATGGTATTTATATC','3','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                    ['sseq4','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
                expected = ["\t".join(x) for x in expected]+['']
                observed = subprocess.check_output(cmd, shell=True).split("\n")
                self.assertEqual(expected, observed)


if __name__ == "__main__":
    unittest.main()
