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

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.appraiser import Appraiser
from singlem.otu_table_collection import OtuTableCollection

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None

    def test_hello_world(self):
        metagenome_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual(7, a.num_found)
        self.assertEqual(4, a.num_not_found)
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(1, len(a.found_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.found_otus[0].sequence)
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                         a.not_found_otus[0].sequence)
        
    def test_multiple_samples(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(2, len(app.appraisal_results))
        a = app.appraisal_results[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(7, a.num_found)
        self.assertEqual(0, a.num_not_found)
        a = app.appraisal_results[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual(0, a.num_found)
        self.assertEqual(4, a.num_not_found)
        
    def test_clusterer_all_cluster_all_good(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATT','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.5)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(11, a.num_found)
        self.assertEqual(0, a.num_not_found)
        
        
    def test_clusterer_all_cluster_some_good(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','minimal',     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.7)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(7, a.num_found)
        self.assertEqual(4, a.num_not_found)
        
        
    def test_clusterer_all_cluster_two_samples(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','maximal',     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.7)
        self.assertEqual(2, len(app.appraisal_results))
        a = app.appraisal_results[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(7, a.num_found)
        self.assertEqual(0, a.num_not_found)
        self.assertEqual(1, len(a.found_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.found_otus[0].sequence)
        self.assertEqual(0, len(a.not_found_otus))
        a = app.appraisal_results[0]
        self.assertEqual('maximal', a.metagenome_sample_name)
        self.assertEqual(0, a.num_found)
        self.assertEqual(4, a.num_not_found)
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                         a.not_found_otus[0].sequence)
        self.assertEqual(0, len(a.found_otus))
        
        
    def test_clusterer_all_cluster_two_samples_some_cluster(self):
        # non-As and genome cluster together but are not exactly the same
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.12.ribosomal_protein_L11_rplK','minimal','AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','12','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','maximal',     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
                    ['4.11.ribosomal_protein_L10','maximal',     'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATG','1','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.7)
        self.assertEqual(2, len(app.appraisal_results))
        
        a = app.appraisal_results[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(7, a.num_found)
        self.assertEqual(12, a.num_not_found)
        self.assertEqual(1, len(a.found_otus))
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.found_otus[0].sequence)
        self.assertEqual('minimal',
                         a.found_otus[0].sample_name)
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                         a.not_found_otus[0].sequence)
        self.assertEqual('minimal',
                         a.not_found_otus[0].sample_name)
        
        a = app.appraisal_results[0]
        self.assertEqual('maximal', a.metagenome_sample_name)
        self.assertEqual(1, a.num_found)
        self.assertEqual(4, a.num_not_found)
        self.assertEqual(1, len(a.found_otus))
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATG',
                         a.found_otus[0].sequence)
        self.assertEqual('maximal',
                         a.found_otus[0].sample_name)
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                         a.not_found_otus[0].sequence)
        self.assertEqual('maximal',
                         a.not_found_otus[0].sample_name)

    def test_print_appraisal(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        
        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(2, len(app.appraisal_results))
        a = app.appraisal_results[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(7, a.num_found)
        self.assertEqual(0, a.num_not_found)
        a = app.appraisal_results[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual(0, a.num_found)
        self.assertEqual(4, a.num_not_found)
        
        to_print = StringIO()
        appraiser.print_appraisal(app, to_print)
        self.assertEqual("sample\tnum_found\tnum_not_found\tpercent_found\nanother\t0\t4\t0.0\nminimal\t7\t0\t100.0\ntotal\t7\t4\t63.6\naverage\t3.5\t2.0\t50.0\n", to_print.getvalue())
        
        to_print = StringIO()
        found_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, to_print,
                                  accounted_for_otu_table_io=found_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io)
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'])
                          ])+"\n",
                         found_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'])])+"\n",
                         not_found_otu_table_io.getvalue())
        
    def test_contamination(self):
        metagenome_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'minimal',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            [
                '4.16.ribosomal_protein_S5',
                'minimal',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '8','12.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        genomes_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC', #one base pair different to the one above
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual(8, a.num_found)
        self.assertEqual(7, a.num_not_found)

                
    def test_contamination_near_enough(self):
        metagenome_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'minimal',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            [
                '4.16.ribosomal_protein_S5',
                'another',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '8','12.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        genomes_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC', #one base pair different to the one above
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGGGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',#one base pair different to the one above
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection)
        self.assertEqual(2, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual(0, a.num_found)
        self.assertEqual(8, a.num_not_found)
        a = app.appraisal_results[1]
        self.assertEqual(0, a.num_found)
        self.assertEqual(7, a.num_not_found)

        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.9)
        self.assertEqual(2, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual(8, a.num_found)
        self.assertEqual(0, a.num_not_found)
        a = app.appraisal_results[1]
        self.assertEqual(0, a.num_found)
        self.assertEqual(7, a.num_not_found)


       

if __name__ == "__main__":
    unittest.main()
