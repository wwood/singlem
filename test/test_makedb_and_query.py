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
import itertools

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    query_result_headers = ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy']
    protein_query_result_headers = ['query_name','query_sequence','query_protein_sequence','divergence','num_hits','sample','marker','hit_sequence','hit_protein_sequence','taxonomy']
    maxDiff = None

    def assertEqualOtuTable(self, expected_array, observed_string, message=None):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0], message)

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]), message)

    def test_makedb_and_dump(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = "%s makedb --db_path %s/db --otu_table %s/methanobacteria/otus.transcripts.on_target.csv --sequence-database-methods none" %(path_to_script,
                                                            d,
                                                            path_to_data)
            extern.run(cmd)

            observed = extern.run("%s query --dump --db %s/db" % (path_to_script, d))
            with open('%s/methanobacteria/otus.transcripts.on_target.csv' % path_to_data) as f:
                expected = [l.split("\t") for l in f.read().split("\n")]
            self.assertEqualOtuTable(expected, observed)

    def test_makedb_query_methanobacteria(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = "%s makedb --db %s/db --otu-table %s/methanobacteria/otus.transcripts.on_target.csv --sequence-database-methods annoy scann nmslib naive" %(
                path_to_script,
                d,
                path_to_data)
            extern.run(cmd)

            for method in ['annoy','nmslib','scann','naive']:
                cmd = "%s query --query-otu-table %s/methanobacteria/otus.transcripts.on_target.3random.csv --db %s/db --search-method %s --max-nearest-neighbours 2" % (
                    path_to_script,
                    path_to_data,
                    d,
                    method)
                observed = extern.run(cmd)
                self.assertEqual(observed.split("\n")[0], "\t".join(self.query_result_headers))
                self.assertTrue('GB_GCA_000309865.1_protein	CAGACTGAAATATTCATGGACAACATGCGAATGTTCCTTAAAGAAGAGGGCCAGGGGATG	0	1	GB_GCA_000309865.1_protein	S3.32.Fibrillarin	CAGACTGAAATATTCATGGACAACATGCGAATGTTCCTTAAAGAAGAGGGCCAGGGGATG	Root; d__Archaea; p__Methanobacteriota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__Methanobacterium sp000309865\n' in observed)

    def test_protein_search_methanobacteria(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = "%s makedb --db %s/db --otu-table %s/methanobacteria/otus.transcripts.on_target.csv --sequence-database-methods annoy scann nmslib naive" %(path_to_script,
                                                            d,
                                                            path_to_data)
            extern.run(cmd)

            for method in ['annoy','nmslib','scann','naive']:
                cmd = "%s query --sequence-type protein --query-otu-table %s/methanobacteria/otus.transcripts.on_target.3random.csv --db %s/db --search-method %s --max-nearest-neighbours 2" % (
                    path_to_script,
                    path_to_data,
                    d,
                    method)
                observed = extern.run(cmd)
                self.assertEqual(observed.split("\n")[0], "\t".join(self.protein_query_result_headers))
                self.assertTrue('GB_GCA_000309865.1_protein	CAGACTGAAATATTCATGGACAACATGCGAATGTTCCTTAAAGAAGAGGGCCAGGGGATG	QTEIFMDNMRMFLKEEGQGM	0	1	RS_GCF_000302455.1_protein	S3.32.Fibrillarin	CAGACTGAAATATTCATGGACAACATGCGAATGTTCCTGAAAGAAGAGGGTCAGGGAATG	QTEIFMDNMRMFLKEEGQGM	Root; d__Archaea; p__Methanobacteriota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__Methanobacterium formicicum_A\n' in observed)


    def test_protein_search_methanobacteria_preload_db(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = "%s makedb --db %s/db --otu-table %s/methanobacteria/otus.transcripts.on_target.csv --sequence-database-methods annoy scann nmslib naive" %(path_to_script,
                                                            d,
                                                            path_to_data)
            extern.run(cmd)

            for method in ['scann','naive']: # not implemented for 'annoy','nmslib',
                cmd = "%s query --preload-db --sequence-type protein --query-otu-table %s/methanobacteria/otus.transcripts.on_target.3random.csv --db %s/db --search-method %s --max-nearest-neighbours 2" % (
                    path_to_script,
                    path_to_data,
                    d,
                    method)
                observed = extern.run(cmd)
                self.assertEqual(observed.split("\n")[0], "\t".join(self.protein_query_result_headers))
                self.assertTrue('GB_GCA_000309865.1_protein	CAGACTGAAATATTCATGGACAACATGCGAATGTTCCTTAAAGAAGAGGGCCAGGGGATG	QTEIFMDNMRMFLKEEGQGM	0	1	RS_GCF_000302455.1_protein	S3.32.Fibrillarin	CAGACTGAAATATTCATGGACAACATGCGAATGTTCCTGAAAGAAGAGGGTCAGGGAATG	QTEIFMDNMRMFLKEEGQGM	Root; d__Archaea; p__Methanobacteriota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__Methanobacterium formicicum_A\n' in observed)





    def test_limit_per_sequence(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = "%s makedb --db_path %s/db --otu_table %s/methanobacteria/otus.transcripts.on_target.csv --sequence-database-methods annoy scann nmslib naive" %(path_to_script,
                                                            d,
                                                            path_to_data)
            extern.run(cmd)

            for (seq_type,unlimited_count,limited_count) in [('nucleotide',64,61),('protein',72,39)]:
                for method in ['annoy','nmslib','naive','scann']:

                    cmd = "%s query --sequence-type %s --query-otu-table %s/methanobacteria/otus.transcripts.on_target.3random.csv --db %s/db --search-method %s" % (
                        path_to_script,
                        seq_type,
                        path_to_data,
                        d,
                        method)
                    observed = extern.run(cmd)
                    self.assertEqual(unlimited_count+1, len(observed.split("\n")), 'umlimited %s %s' %(seq_type,method))

                    cmd = "%s query --limit-per-sequence 1 --sequence-type %s --query-otu-table %s/methanobacteria/otus.transcripts.on_target.3random.csv --db %s/db --search-method %s" % (
                        path_to_script,
                        seq_type,
                        path_to_data,
                        d,
                        method)
                    observed = extern.run(cmd)
                    self.assertEqual(limited_count+1, len(observed.split("\n")), 'limited %s %s' %(seq_type,method))

    def test_query_with_otu_table_two_samples(self):
        with tempfile.NamedTemporaryFile(mode='w') as f:
            query = [self.headers,
                     # second sequence with an extra A at the end
                     ['ribosomal_protein_L11_rplK_gpkg','maximal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
                     ] # converted A to T in the middle
            query = "\n".join(["\t".join(x) for x in query])
            f.write(query)
            f.flush()

            cmd = "%s query --query-otu-table %s --db %s" % (
                path_to_script,
                f.name,
                os.path.join(path_to_data,'a.sdb'))

            expected = 'query_name\tquery_sequence\tdivergence\tnum_hits\tsample\tmarker\thit_sequence\ttaxonomy\nminimal\tCGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA\t40\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\nmaximal\tCGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA\t40\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\n'.split('\n')
            observed = extern.run(cmd).split('\n')
            self.assertEqual(expected, observed)



    def test_query_with_otu_table_two_samples_preload_db(self):
        with tempfile.NamedTemporaryFile(mode='w') as f:
            query = [self.headers,
                     # second sequence with an extra A at the end
                     ['ribosomal_protein_L11_rplK_gpkg','maximal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                     ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
                     ] # converted A to T in the middle
            query = "\n".join(["\t".join(x) for x in query])
            f.write(query)
            f.flush()

            cmd = "%s query --preload-db --query-otu-table %s --db %s" % (
                path_to_script,
                f.name,
                os.path.join(path_to_data,'a.sdb'))

            expected = 'query_name\tquery_sequence\tdivergence\tnum_hits\tsample\tmarker\thit_sequence\ttaxonomy\nminimal\tCGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA\t40\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\nmaximal\tCGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA\t40\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\n'.split('\n')
            observed = extern.run(cmd).split('\n')
            self.assertEqual(expected, observed)



    # def test_query_with_otu_table_two_samples_same_sequence(self):
    #     with tempfile.NamedTemporaryFile(mode='w') as f:
    #         query = [self.headers,
    #                  # second sequence with an extra A at the end
    #                  ['ribosomal_protein_L11_rplK_gpkg','maximal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                  ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli']
    #                  ] # converted A to T in the middle
    #         query = "\n".join(["\t".join(x) for x in query])
    #         f.write(query)
    #         f.flush()

    #         with tempfile.TemporaryDirectory() as d:
    #             cmd = "{} makedb --db {}/sdb --otu_table {}".format(
    #                 path_to_script, d, f.name)
    #             extern.run(cmd)

    #             cmd = "{} query --query_otu_table {} --db {}/sdb".format(
    #                 path_to_script,
    #                 f.name,
    #                 d)

    #             expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
    #                         ['maximal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','maximal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                         ['maximal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #                         ['minimal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','maximal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                         ['minimal;ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #                         ]
    #             observed = extern.run(cmd)
    #             self.assertEqualOtuTable(expected, observed)

    # def test_n_vs_gap(self):
    #     otu_table = [self.headers,
    #                  # same as query
    #                  ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #                  # Ns instead of gaps, plus and A instead of C at end
    #                  ['ribosomal_protein_L11_rplK_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAANNNTATATCTTCACTGAGAGAAATGGTATTTATATA','7','4.95','Root; k__Bacteria; p__Firmicutes']]
    #     otu_table = "\n".join(["\t".join(x) for x in otu_table])

    #     with tempfile.NamedTemporaryFile(mode='w') as f:
    #         f.write(otu_table)
    #         f.flush()

    #         with tempfile.TemporaryDirectory() as d:
    #             #d = '/tmp/d'
    #             cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
    #                                                             d,
    #                                                             f.name)
    #             subprocess.check_call(cmd, shell=True)

    #             cmd = "%s query --query_sequence %s --db %s/db" % (
    #                 path_to_script,
    #                 'CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC',
    #                 d)

    #             expected = [['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
    #                         ['unnamed_sequence','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #                         ['unnamed_sequence','CGTCGTTGGAACCCAAAAATGAAA---TATATCTTCACTGAGAGAAATGGTATTTATATC','4','7','minimal','ribosomal_protein_L11_rplK_gpkg','CGTCGTTGGAACCCAAAAATGAAANNNTATATCTTCACTGAGAGAAATGGTATTTATATA','Root; k__Bacteria; p__Firmicutes']]
    #             expected = ["\t".join(x) for x in expected]+['']
    #             self.assertEqual(expected,
    #                              extern.run(cmd).split('\n'))

    # def test_duplicate_seqs_then_another(self):
    #     '''This tests when a two samples have the same OTU, and then there's another separate OTU after that'''
    #     otu_table = [self.headers,
    #                  ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                  ['ribosomal_protein_L11_rplK_gpkg','maximal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                  ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #                  ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
    #     otu_table = "\n".join(["\t".join(x) for x in otu_table])

    #     with tempfile.NamedTemporaryFile(mode='w') as f:
    #         f.write(otu_table)
    #         f.flush()

    #         with tempfile.TemporaryDirectory() as d:
    #             #d = '/tmp/d'
    #             cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
    #                                                             d,
    #                                                             f.name)
    #             subprocess.check_call(cmd, shell=True)

    #             with tempfile.NamedTemporaryFile(mode='w') as infasta:
    #                 infasta.write(">1_\n") # same as the first sequence
    #                 infasta.write("GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\n")
    #                 infasta.write(">2_\n") # third sequence with an replacement A at the end
    #                 infasta.write("CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA\n")
    #                 infasta.flush()

    #                 cmd = "%s query --query_fasta %s --db %s/db" % (
    #                     path_to_script,
    #                     infasta.name,
    #                     d)

    #                 expected = [
    #                     ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
    #                     ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                     ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','9','maximal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                     ['2_','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATA','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
    #                 expected = ["\t".join(x) for x in expected]+['']
    #                 self.assertEqual(expected,
    #                                  extern.run(cmd).split("\n"))


    # def test_divergence0(self):
    #     '''Zero divergence is special because it uses the sqlite DB to query, not blast'''
    #     otu_table = [self.headers,
    #                  ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                  ['ribosomal_protein_L11_rplK_gpkg','maximal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','9','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                  ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','4.95','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
    #     otu_table = "\n".join(["\t".join(x) for x in otu_table])



    #     with tempfile.NamedTemporaryFile(mode='w') as f:
    #         f.write(otu_table)
    #         f.flush()

    #         with tempfile.TemporaryDirectory() as d:

    #             cmd = "%s makedb --db_path %s/db --otu_table %s" %(path_to_script,
    #                                                             d,
    #                                                             f.name)
    #             subprocess.check_call(cmd, shell=True)

    #             with tempfile.NamedTemporaryFile(mode='w') as infasta:
    #                 infasta.write(">1_\n") # same as the first sequence
    #                 infasta.write("GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\n")
    #                 infasta.write(">2_\n") # same as first
    #                 infasta.write("GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\n")
    #                 infasta.write(">3_\n") # same as third in OTU table
    #                 infasta.write("CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC\n")
    #                 infasta.flush()

    #                 cmd = "%s query --query-otu-table %s --db %s/db --max_divergence 0" % (
    #                     path_to_script,
    #                     f.name,
    #                     d)

    #                 expected = [
    #                     ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
    #                     ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                     ['1_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','9','maximal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                     ['2_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','7','minimal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                     ['2_','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','0','9','maximal','ribosomal_protein_L11_rplK_gpkg','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #                     ['3_','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','0','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
    #                 expected = ["\t".join(x) for x in expected]+['']
    #                 self.assertEqual(sorted(expected),
    #                                  sorted(extern.run(cmd).split("\n")))

    # def test_query_by_sample(self):
    #     expected = [
    #         self.headers,
    #         ['ribosomal_protein_L11_rplK_gpkg','m2','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.10','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #         ['ribosomal_protein_S2_rpsB_gpkg','m2','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.40','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #         ['ribosomal_protein_S17_gpkg','m2','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGA','9','20.50','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
    #         ['ribosomal_protein_L11_rplK_gpkg','m3','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.10','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #         ['ribosomal_protein_S2_rpsB_gpkg','m3','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.40','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #         ['ribosomal_protein_S17_gpkg','m3','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','19.50','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
    #     expected = ["\t".join(x) for x in expected]+['']

    #     cmd = "%s query --db %s/b.sdb --sample_names m2 m3" %(path_to_script,
    #                                                           path_to_data)
    #     self.assertEqual(expected, extern.run(cmd).split('\n'))

    # def test_query_by_taxonomy(self):
    #     expected = [
    #         self.headers,
    #         ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.10','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #         ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','19.50','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
    #     expected = ["\t".join(x) for x in expected]+['']

    #     cmd = "%s query --db %s/a.sdb --taxonomy o__Bacillales" %(path_to_script,
    #                                                               path_to_data)
    #     self.assertEqual(expected, extern.run(cmd).split('\n'))

    # def test_query_subject_otu_tables(self):
    #     with tempfile.NamedTemporaryFile(mode='w') as f:
    #         query = "\n".join([">seq1 comment",'CGTCGTTGGAACCCAAAAATGAAAAAATATATCTatgTCACTGAGAGAAATGGTATTTATATC',
    #                            ">sseq4",       'CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC',''])
    #         f.write(query)
    #         f.flush()

    #         subject = [
    #             self.headers,
    #             ['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','15.1','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
    #             ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','12.4','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #             ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','19.5','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
    #         with tempfile.NamedTemporaryFile(mode='w') as subject_f:
    #             subject_f.write("\n".join(["\t".join(x) for x in subject]+['']))
    #             subject_f.flush()

    #             cmd = "%s query --query_fasta %s --subject_otu_tables %s" % (
    #                 path_to_script,
    #                 f.name,
    #                 subject_f.name)

    #             expected = [
    #                 ['query_name','query_sequence','divergence','num_hits','sample','marker','hit_sequence','taxonomy'],
    #                 ['seq1','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTatgTCACTGAGAGAAATGGTATTTATATC','3','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
    #                 ['sseq4','CGTCGTTGGAACCCAAAAATGAAATAATATATCTTCACTGAGAGAAATGGTATTTATATC','1','6','minimal','ribosomal_protein_S2_rpsB_gpkg','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','Root; k__Bacteria; p__Firmicutes; c__Bacilli']]
    #             expected = ["\t".join(x) for x in expected]+['']
    #             observed = extern.run(cmd).split("\n")
    #             self.assertEqual(expected, observed)

if __name__ == "__main__":
    unittest.main()
