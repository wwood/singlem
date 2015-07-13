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
import sys
import os
from StringIO import StringIO

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.query_formatters import SparseResultFormatter, DenseResultFormatter,\
    NamedQueryDefinition, NameSequenceQueryDefinition
from singlem.sequence_database import DBSequence

class Tests(unittest.TestCase):
    subjects = []
    
    def setUp(self):
        self.s = []
        for i, row in enumerate([['ribosomal_protein_L11_rplK_gpkg','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['ribosomal_protein_S2_rpsB_gpkg','minimal','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','6','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
                    ['ribosomal_protein_S17_gpkg','minimal','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','9','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]):
            sub = DBSequence()
            sub.sequence_id = i
            sub.marker = row[0]
            sub.sample_name = row[1]
            sub.sequence = row[2]
            sub.count = int(row[3])
            sub.taxonomy = row[4]
            self.s.append(sub)

    def test_sparse(self):
        div_subs = [[2,self.s[0]],
                    [1,self.s[1]]] 
        formatter = SparseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(['query_name\tdivergence\tnum_hits\tsample\tmarker\thit_sequence\ttaxonomy',
         'sampleme1\t1\t6\tminimal\tribosomal_protein_S2_rpsB_gpkg\tCGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli',
         'sampleme1\t2\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales',
         ''], strio.getvalue().split("\n"))
        
    def test_sparse_query_sequence_definition(self):
        div_subs = [[2,self.s[0]],
                    [1,self.s[1]]] 
        formatter = SparseResultFormatter(NameSequenceQueryDefinition(['sampleme1'],['ATGC']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(['query_name\tquery_sequence\tdivergence\tnum_hits\tsample\tmarker\thit_sequence\ttaxonomy',
         'sampleme1\tATGC\t1\t6\tminimal\tribosomal_protein_S2_rpsB_gpkg\tCGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli',
         'sampleme1\tATGC\t2\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales',
         ''], strio.getvalue().split("\n"))
        
    def test_sparse_two_samples(self):
        div_subs = [[[2,self.s[0]],
                     [1,self.s[1]]],
                    [[4,self.s[1]],
                     [5,self.s[2]]]] 
        formatter = SparseResultFormatter(NamedQueryDefinition(['sampleme1','sam2']), div_subs)
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(['query_name\tdivergence\tnum_hits\tsample\tmarker\thit_sequence\ttaxonomy',
         'sampleme1\t1\t6\tminimal\tribosomal_protein_S2_rpsB_gpkg\tCGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli',
         'sampleme1\t2\t7\tminimal\tribosomal_protein_L11_rplK_gpkg\tGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales',
         'sam2\t4\t6\tminimal\tribosomal_protein_S2_rpsB_gpkg\tCGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli',
         'sam2\t5\t9\tminimal\tribosomal_protein_S17_gpkg\tGCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT\tRoot; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus',
         ''], strio.getvalue().split("\n"))
        
    def test_dense_single_sample(self):
        div_subs = [[2,self.s[0]]] 
        formatter = DenseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(["\t".join(x) for x in [['query_name','hit_sequence','divergence','minimal','taxonomy'],
                          ['sampleme1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','2','7','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],['']]],
                         strio.getvalue().split("\n"))
        
    def test_dense_two_seqs_same_sample(self):
        div_subs = [[2,self.s[0]],
                    [1,self.s[1]]] 
        formatter = DenseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(["\t".join(x) for x in [['query_name','hit_sequence','divergence','minimal','taxonomy'],
             ['sampleme1','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','1','6','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
             ['sampleme1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','2','7','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],['']]],
                         strio.getvalue().split("\n"))
        
    def test_dense_two_seqs_diff_sample(self):
        self.s[1].sample_name='m2'
        div_subs = [[2,self.s[0]],
                    [1,self.s[1]]] 
        formatter = DenseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(["\t".join(x) for x in [['query_name','hit_sequence','divergence','m2','minimal','taxonomy'],
             ['sampleme1','CGTCGTTGGAACCCAAAAATGAAAAAATATATCTTCACTGAGAGAAATGGTATTTATATC','1','6','0','Root; k__Bacteria; p__Firmicutes; c__Bacilli'],
             ['sampleme1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','2','0','7','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],['']]],
                         strio.getvalue().split("\n"))
        
    def test_dense_one_seq_two_samples(self):
        self.s[1].sequence = self.s[0].sequence
        self.s[1].sample_name='m2'
        div_subs = [[1,self.s[0]],
                    [1,self.s[1]]] 
        formatter = DenseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(["\t".join(x) for x in [['query_name','hit_sequence','divergence','minimal','m2','taxonomy'],
             ['sampleme1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','7','6','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],['']]],
                         strio.getvalue().split("\n"))
        
    def test_dense_two_seq_three_samples_order1(self):
        self.s[1].sequence = self.s[0].sequence
        self.s[1].sample_name='m2'
        self.s[2].sample_name='m3'
        div_subs = [[1,self.s[0]],
                    [1,self.s[1]],
                    [1,self.s[2]]] 
        formatter = DenseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(["\t".join(x) for x in [['query_name','hit_sequence','divergence','minimal','m2','m3','taxonomy'],
             ['sampleme1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','7','6','0','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['sampleme1','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','1','0','0','9','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
             ['']]],strio.getvalue().split("\n"))
        
    def test_dense_two_seq_three_samples_order2(self):
        self.s[1].sequence = self.s[0].sequence
        self.s[1].sample_name='m2'
        self.s[2].sample_name='m3'
        self.s[2].count = 100
        div_subs = [[1,self.s[0]],
                    [1,self.s[1]],
                    [1,self.s[2]]] 
        formatter = DenseResultFormatter(NamedQueryDefinition(['sampleme1']), [div_subs])
        strio = StringIO()
        formatter.write(strio)
        self.assertEqual(["\t".join(x) for x in [['query_name','hit_sequence','divergence','m3','minimal','m2','taxonomy'],
             ['sampleme1','GCTAAATTAGGAGACATTGTTAAAATTCAAGAAACTCGTCCTTTATCAGCAACAAAACGT','1','100','0','0','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
             ['sampleme1','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','0','7','6','Root; k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
             ['']]],strio.getvalue().split("\n"))
                            
if __name__ == "__main__":
    unittest.main()
