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
import sys
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.otu_table import OtuTable
from singlem.condense  import Condenser
from singlem.archive_otu_table import ArchiveOtuTable
from singlem.pipe import QUERY_BASED_ASSIGNMENT_METHOD, DIAMOND_ASSIGNMENT_METHOD

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data','condense')

class Tests(unittest.TestCase):
    maxDiff = None

    def test_apply_expectation_maximization_core_trivial(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 1}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.05,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        print(species_to_coverage)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 1.05},
            species_to_coverage
        )
    
    def test_apply_expectation_maximization_core_split1(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 1, 'Root;d__Viruses;p;c;o;f;g;tax2': 1}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1','Root;d__Viruses;p;c;o;f;g;tax2'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split2(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 1, 'Root;d__Viruses;p;c;o;f;g;tax2': 1}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1','Root;d__Viruses;p;c;o;f;g;tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 2.2,},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split3(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 1, 'Root;d__Viruses;p;c;o;f;g;tax2': 1, 'Root;d__Viruses;p;c;o;f;g;tax3': 1}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1','Root;d__Viruses;p;c;o;f;g;tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1','Root;d__Viruses;p;c;o;f;g;tax2','Root;d__Viruses;p;c;o;f;g;tax3'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 1.1, 'Root;d__Viruses;p;c;o;f;g;tax2': 1.1},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split4(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 1, 'Root;d__Viruses;p;c;o;f;g;tax2': 1, 'Root;d__Viruses;p;c;o;f;g;tax3': 1, 'Root;d__Viruses;p;c;o;f;g;tax4': 1, 'Root;d__Viruses;p;c;o;f;g;tax5': 1}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1','Root;d__Viruses;p;c;o;f;g;tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1','Root;d__Viruses;p;c;o;f;g;tax2','Root;d__Viruses;p;c;o;f;g;tax3'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.2,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax5','Root;d__Viruses;p;c;o;f;g;tax4'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 1.1, 'Root;d__Viruses;p;c;o;f;g;tax2': 1.1, 'Root;d__Viruses;p;c;o;f;g;tax4': 0.6, 'Root;d__Viruses;p;c;o;f;g;tax5': 0.6},
            species_to_coverage
        )
        self.assertEqual(sorted([
            ['Root;d__Viruses;p;c;o;f;g;tax1', 'Root;d__Viruses;p;c;o;f;g;tax2'],
            ['Root;d__Viruses;p;c;o;f;g;tax1', 'Root;d__Viruses;p;c;o;f;g;tax2', 'Root;d__Viruses;p;c;o;f;g;tax3'],
            ['Root;d__Viruses;p;c;o;f;g;tax5', 'Root;d__Viruses;p;c;o;f;g;tax4']
        ]), sorted(best_hit_taxonomy_sets))

    def test_apply_expectation_maximization_core_trimmed_mean_trivial(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 3}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.06,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g2 sample1 seq1')+[1,1.05,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0.4, {'Viruses': ['g1','g2','g3']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_trimmed_mean_two_species(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 3, 'Root;d__Viruses;p;c;o;f;g;tax2': 3}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.06,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g2 sample1 seq1')+[1,1.05,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1', 'Root;d__Viruses;p;c;o;f;g;tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0.4, {'Viruses': ['g1','g2','g3']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;g;tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_trimmed_mean_no_hits(self):
        taxon_marker_counts = {'Root;d__Viruses;p;c;o;f;g;tax1': 3, 'Root;d__Viruses;p;c;o;f;g;tax2': 3}
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            # same gene for each OTU
            str.split('g1 sample1 seq1')+[1,1.06,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.05,'','','','','',['Root;d__Viruses;p;c;o;f;g;tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0.4, {'Viruses': ['g1','g2','g3']}, min_genes_for_whitelist=0, proximity_cutoff=0, taxon_marker_counts=taxon_marker_counts)
        self.assertEqual(
            {},
            species_to_coverage
        )
    
    def test_apply_genus_expectation_maximization_core_trivial(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.05,'','','','','',['Root;d__Viruses;p;c;o;f;tax1'],DIAMOND_ASSIGNMENT_METHOD]
        ]
        best_hit_taxonomy_sets = Condenser()._apply_genus_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, avg_num_genes_per_species=1)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;tax1': 1.05},
            best_hit_taxonomy_sets[0]
        )
    
    def test_apply_genus_expectation_maximization_core_split1(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;tax1'],DIAMOND_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root;d__Viruses;p;c;o;f;tax1','Root;d__Viruses;p;c;o;f;tax2'],DIAMOND_ASSIGNMENT_METHOD]
        ]
        best_hit_taxonomy_sets = Condenser()._apply_genus_expectation_maximization_core(otus, 0, {'Viruses': ['g1']}, avg_num_genes_per_species=3)
        self.assertEqual(
            {'Root;d__Viruses;p;c;o;f;tax1': 0.733,
             'Root;d__Viruses;p;c;o;f;tax2': 0.001},
            best_hit_taxonomy_sets[0]
        )

if __name__ == "__main__":
    import logging
    # logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    unittest.main()