#!/usr/bin/env python3

#=======================================================================
# Authors: Rossen Zhao
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
import tempdir
import sys
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.otu_table import OtuTable
from singlem.condense  import Condenser

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data','condense')

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_condense_small(self): #TODO add singlem packages, example input and output comparators
        with tempdir.in_tempdir():
            # trim 35% so 1 trimmed off top and bottom of 3 spkgs
            cmd = "{} condense --input-otu-tables {}/small_condense_input.csv --trim-percent 35 --output-otu-table small_condense_output.csv --singlem-packages {}/*spkg".format(
                path_to_script, path_to_data, path_to_data, path_to_data
            )
            observed = extern.run(cmd)

            with open('small_condense_output.csv') as observed:
                with open(os.path.join(path_to_data, 'small_condense_output.csv')) as expected:
                    self.assertListEqual(list(expected), list(observed))
    
    def test_condense_zero_trim(self): #TODO add singlem packages, example input and output comparators
        with tempdir.in_tempdir():
            cmd = "{} condense --input-otu-tables {}/small_condense_input.csv --trim-percent 0 --output-otu-table small_condense_output.csv --singlem-packages {}/*spkg".format(
                path_to_script, path_to_data, path_to_data, path_to_data
            )
            observed = extern.run(cmd)

            with open('small_condense_output.csv') as observed:
                with open(os.path.join(path_to_data, 'small_condense_output_no_trim.csv')) as expected:
                    self.assertListEqual(list(expected), list(observed))
    
    def test_condense_two_tables_cmdline(self): #TODO add singlem packages, example input and output comparators
        with tempdir.in_tempdir():
            cmd = "{} condense --input-otu-tables {}/small_condense_input_another.csv --trim-percent 35 --output-otu-table out --krona a.html --singlem-packages {}/*spkg".format(
                path_to_script, path_to_data, path_to_data, path_to_data
            )
            observed = extern.run(cmd)
            
            with open('out') as observed:
                with open(os.path.join(path_to_data, 'small_condense_output_2_samples.csv')) as expected:
                    self.assertListEqual(list(observed), list(expected))

    def test_apply_expectation_maximization_core_trivial(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.05,'tax1',['tax1']]
        ]
        species_to_coverage, best_hit_taxonomy_sets, best_hits_field_index = Condenser()._apply_expectation_maximization_core(otus)
        self.assertEqual(
            {'tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split1(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1']],
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']]
        ]
        species_to_coverage, best_hit_taxonomy_sets, best_hits_field_index = Condenser()._apply_expectation_maximization_core(otus)
        self.assertEqual(
            {'tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split2(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1']]
        ]
        species_to_coverage, best_hit_taxonomy_sets, best_hits_field_index = Condenser()._apply_expectation_maximization_core(otus)
        self.assertEqual(
            {'tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split3(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        species_to_coverage, best_hit_taxonomy_sets, best_hits_field_index = Condenser()._apply_expectation_maximization_core(otus)
        self.assertEqual(
            {'tax1': 1.1, 'tax2': 1.1},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split4(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']],
            str.split('g1 sample1 seq1')+[1,1.2,'',['tax5','tax4']]
        ]
        species_to_coverage, best_hit_taxonomy_sets, best_hits_field_index = Condenser()._apply_expectation_maximization_core(otus)
        self.assertEqual(
            {'tax1': 1.1, 'tax2': 1.1, 'tax4': 0.6, 'tax5': 0.6},
            species_to_coverage
        )
        self.assertEqual(sorted([
            ['tax1', 'tax2'],
            ['tax1', 'tax2', 'tax3'],
            ['tax5', 'tax4']
        ]), sorted(best_hit_taxonomy_sets))
        self.assertEqual(6, best_hits_field_index)

    def test_gather_equivalence_classes_from_list_of_species_lists1(self):
        species_lists = [['tax1'], ['tax2']]
        expected = {
            'tax1': {'tax1'},
            'tax2': {'tax2'}
        }
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_species_lists(species_lists)
        )

    def test_gather_equivalence_classes_from_list_of_species_lists2(self):
        species_lists = [['tax1'], ['tax2'],['tax1','tax2']]
        expected = {
            'tax1': {'tax1'},
            'tax2': {'tax2'}
        }
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_species_lists(species_lists)
        )

    def test_gather_equivalence_classes_from_list_of_species_lists3(self):
        species_lists = [['tax1','tax2'],['tax1'], ['tax2']]
        expected = {
            'tax1': {'tax1'},
            'tax2': {'tax2'}
        }
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_species_lists(species_lists)
        )

    def test_gather_equivalence_classes_from_list_of_species_lists4(self):
        species_lists = [['tax1','tax2'],['tax1','tax2','tax3']]
        expected = {'tax1': {'tax1', 'tax2'}, 'tax2': {'tax1', 'tax2'}, 'tax3': {'tax3'}}
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_species_lists(species_lists)
        )

    def test_demultiplex_otus1(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['tax1']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        expected = [
            str.split('g1 sample1 seq1')+[1,1.1,'tax1'],
            # str.split('g1 sample1 seq1')+[1,1.1,'tax2']
        ]
        self.assertEqual(
            expected,
            Condenser()._demultiplex_otus(otus, {'tax1': 1}, 6, {'tax1': {'tax1'}}).data
        )

    def test_demultiplex_otus2(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['Rooter; tax1','Rooter; tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        expected = [
            str.split('g1 sample1 seq1')+[1,1.1,'Rooter'],
            # str.split('g1 sample1 seq1')+[1,1.1,'tax2']
        ]
        self.assertEqual(
            expected,
            Condenser()._demultiplex_otus(otus, {'Rooter; tax1': 1, 'Rooter; tax2': 1}, 6, \
                {'Rooter; tax1': {'Rooter; tax1','Rooter; tax2'}, 'Rooter; tax2': {'Rooter; tax1','Rooter; tax2'}} \
                ).data
        )

    def test_demultiplex_otus3(self):
        otus = OtuTable()
        otus.fields = str.split('gene sample sequence num_hits coverage taxonomy equal_best_hit_taxonomies')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'',['Rooter; tax1','Rooter; tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        expected = [
            str.split('g1 sample1 seq1')+[1,0.55,'Rooter; tax1'],
            str.split('g1 sample1 seq1')+[1,0.55,'Rooter; tax2'],
            # str.split('g1 sample1 seq1')+[1,1.1,'tax2']
        ]
        self.assertEqual(
            expected,
            Condenser()._demultiplex_otus(otus, {'Rooter; tax1': 1, 'Rooter; tax2': 1}, 6, \
                {'Rooter; tax1': {'Rooter; tax1'}, 'Rooter; tax2': {'Rooter; tax2'}} \
                ).data
        )

if __name__ == "__main__":
    # import logging
    # logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    unittest.main()