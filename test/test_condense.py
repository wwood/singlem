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
from singlem.condense  import Condenser, WordNode, CondensedCommunityProfile, SylphHit, SylphProfile
from singlem.condense import _canonical_species_key, _gtdb_string_to_wordnode_array
from singlem.archive_otu_table import ArchiveOtuTable
from singlem.pipe import QUERY_BASED_ASSIGNMENT_METHOD

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data','condense')

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_apply_expectation_maximization_core_trivial(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.05,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split1(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split2(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split3(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2','Root; d__Bacteria; p;c;o;f;g; tax3'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.1, 'Root; d__Bacteria; p;c;o;f;g; tax2': 1.1},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split4(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2','Root; d__Bacteria; p;c;o;f;g; tax3'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.2,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax5','Root; d__Bacteria; p;c;o;f;g; tax4'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.1, 'Root; d__Bacteria; p;c;o;f;g; tax2': 1.1, 'Root; d__Bacteria; p;c;o;f;g; tax4': 0.6, 'Root; d__Bacteria; p;c;o;f;g; tax5': 0.6},
            species_to_coverage
        )
        self.assertEqual(sorted([
            ['Root; d__Bacteria; p;c;o;f;g; tax1', 'Root; d__Bacteria; p;c;o;f;g; tax2'],
            ['Root; d__Bacteria; p;c;o;f;g; tax1', 'Root; d__Bacteria; p;c;o;f;g; tax2', 'Root; d__Bacteria; p;c;o;f;g; tax3'],
            ['Root; d__Bacteria; p;c;o;f;g; tax5', 'Root; d__Bacteria; p;c;o;f;g; tax4']
        ]), sorted(best_hit_taxonomy_sets))

    def test_apply_expectation_maximization_core_trimmed_mean_trivial(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.06,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g2 sample1 seq1')+[1,1.05,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0.4, {'Bacteria': ['g1','g2','g3']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_trimmed_mean_two_species(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.06,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g2 sample1 seq1')+[1,1.05,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1', 'Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0.4, {'Bacteria': ['g1','g2','g3']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_trimmed_mean_no_hits(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            # same gene for each OTU
            str.split('g1 sample1 seq1')+[1,1.06,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,1.05,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0.4, {'Bacteria': ['g1','g2','g3']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {},
            species_to_coverage
        )

    def test_gather_equivalence_classes_from_list_of_taxon_lists1(self):
        species_lists = [['tax1'], ['tax2']]
        expected = {
            'tax1': {'tax1'},
            'tax2': {'tax2'}
        }
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_taxon_lists(species_lists)
        )

    def test_gather_equivalence_classes_from_list_of_taxon_lists2(self):
        species_lists = [['tax1'], ['tax2'],['tax1','tax2']]
        expected = {
            'tax1': {'tax1'},
            'tax2': {'tax2'}
        }
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_taxon_lists(species_lists)
        )

    def test_gather_equivalence_classes_from_list_of_taxon_lists3(self):
        species_lists = [['tax1','tax2'],['tax1'], ['tax2']]
        expected = {
            'tax1': {'tax1'},
            'tax2': {'tax2'}
        }
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_taxon_lists(species_lists)
        )

    def test_gather_equivalence_classes_from_list_of_taxon_lists4(self):
        species_lists = [['tax1','tax2'],['tax1','tax2','tax3']]
        expected = {'tax1': {'tax1', 'tax2'}, 'tax2': {'tax1', 'tax2'}, 'tax3': {'tax3'}}
        self.assertEqual(
            expected,
            Condenser()._gather_equivalence_classes_from_list_of_taxon_lists(species_lists)
        )

    def test_demultiplex_otus1(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        expected = [
            str.split('g1 sample1 seq1')+[1,1.1,'tax1','','','','',['tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            # str.split('g1 sample1 seq1')+[1,1.1,'tax2']
        ]
        self.assertEqual(
            expected,
            Condenser()._demultiplex_otus(otus, {'tax1': 1}, {'tax1': {'tax1'}}, QUERY_BASED_ASSIGNMENT_METHOD).data
        )

    def test_demultiplex_otus2(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Rooter; tax1','Rooter; tax2'], QUERY_BASED_ASSIGNMENT_METHOD],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        expected = [
            str.split('g1 sample1 seq1')+[1,1.1,'Rooter','','','','',['Rooter; tax1','Rooter; tax2'], QUERY_BASED_ASSIGNMENT_METHOD],
            # str.split('g1 sample1 seq1')+[1,1.1,'tax2']
        ]
        self.assertEqual(
            expected,
            Condenser()._demultiplex_otus(otus, {'Rooter; tax1': 1, 'Rooter; tax2': 1}, \
                {'Rooter; tax1': {'Rooter; tax1','Rooter; tax2'}, 'Rooter; tax2': {'Rooter; tax1','Rooter; tax2'}}, QUERY_BASED_ASSIGNMENT_METHOD \
                ).data
        )

    def test_demultiplex_otus3(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            str.split('g1 sample1 seq1')+[1,1.1,'','','','','',['Rooter; tax1','Rooter; tax2'], QUERY_BASED_ASSIGNMENT_METHOD],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2']],
            # str.split('g1 sample1 seq1')+[1,1.1,'',['tax1','tax2','tax3']]
        ]
        expected = [
            str.split('g1 sample1 seq1')+[1,0.55,'Rooter; tax1','','','','',['Rooter; tax1','Rooter; tax2'], QUERY_BASED_ASSIGNMENT_METHOD],
            str.split('g1 sample1 seq1')+[1,0.55,'Rooter; tax2','','','','',['Rooter; tax1','Rooter; tax2'], QUERY_BASED_ASSIGNMENT_METHOD],
            # str.split('g1 sample1 seq1')+[1,1.1,'tax2']
        ]
        self.assertEqual(
            expected,
            Condenser()._demultiplex_otus(otus, {'Rooter; tax1': 1, 'Rooter; tax2': 1}, \
                {'Rooter; tax1': {'Rooter; tax1'}, 'Rooter; tax2': {'Rooter; tax2'}}, QUERY_BASED_ASSIGNMENT_METHOD \
                ).data
        )

    # ---- Regime 3: sylph-only species injection ----

    G_ECOLI = 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia'
    S_ECOLI = G_ECOLI + ';s__Escherichia coli'
    S_SHIGELLA = G_ECOLI + ';s__Shigella flexneri'

    def test_canonical_species_key(self):
        self.assertEqual(
            'd__Bacteria;p__P;s__X',
            _canonical_species_key('Root; d__Bacteria; p__P; s__X'))
        self.assertEqual(
            'd__Bacteria;p__P;s__X',
            _canonical_species_key('d__Bacteria;p__P;s__X'))

    def test_gtdb_string_to_wordnode_array(self):
        self.assertEqual(
            ['Root', 'd__Bacteria', 'p__P', 's__X'],
            _gtdb_string_to_wordnode_array('d__Bacteria;p__P;s__X'))
        self.assertEqual(
            ['Root', 'd__Bacteria', 'p__P', 's__X'],
            _gtdb_string_to_wordnode_array('Root; d__Bacteria; p__P; s__X'))

    def test_fit_alpha_enough_anchors(self):
        # singlem = 2 * eff_cov exactly => alpha == 2.0
        singlem = {'a': 10.0, 'b': 20.0, 'c': 30.0}
        sylph = {'a': SylphHit('a', 5.0), 'b': SylphHit('b', 10.0), 'c': SylphHit('c', 15.0)}
        self.assertAlmostEqual(2.0, Condenser()._fit_alpha(singlem, sylph))

    def test_fit_alpha_too_few_anchors(self):
        # Only two species at >= 10x coverage => default to 1.0
        singlem = {'a': 10.0, 'b': 20.0, 'c': 3.0}
        sylph = {'a': SylphHit('a', 5.0), 'b': SylphHit('b', 10.0), 'c': SylphHit('c', 1.5)}
        self.assertEqual(1.0, Condenser()._fit_alpha(singlem, sylph))

    def _find_node(self, profile, word):
        for node in profile.breadth_first_iter():
            if node.word == word:
                return node
        return None

    def test_inject_sylph_only_species_reconciles_with_genus(self):
        root = WordNode(None, 'Root')
        root.add_words(['Root'] + self.S_ECOLI.split(';'), 5.0)
        # Genus-level novel coverage that injection should draw down
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 3.0)
        profile = CondensedCommunityProfile('sample1', root)

        sylph_hits = {
            _canonical_species_key(self.S_ECOLI): SylphHit(self.S_ECOLI, 9.0),  # already present
            _canonical_species_key(self.S_SHIGELLA): SylphHit(self.S_SHIGELLA, 2.0),  # sylph-only
        }
        Condenser()._inject_sylph_only_species(profile, sylph_hits, alpha=1.0)

        # E. coli already in profile, untouched
        self.assertEqual(5.0, self._find_node(profile, 's__Escherichia coli').coverage)
        # Shigella injected at alpha*eff_cov = 2.0
        self.assertEqual(2.0, self._find_node(profile, 's__Shigella flexneri').coverage)
        # Genus novel coverage drawn down from 3.0 to 1.0 (residual 0, no new total)
        self.assertEqual(1.0, self._find_node(profile, 'g__Escherichia').coverage)

    def test_inject_sylph_only_species_residual_is_new_coverage(self):
        root = WordNode(None, 'Root')
        # Only 0.5 of genus novel budget available
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 0.5)
        profile = CondensedCommunityProfile('sample1', root)
        total_before = sum([n.coverage for n in profile.breadth_first_iter()])

        sylph_hits = {_canonical_species_key(self.S_SHIGELLA): SylphHit(self.S_SHIGELLA, 2.0)}
        Condenser()._inject_sylph_only_species(profile, sylph_hits, alpha=1.0)

        self.assertEqual(2.0, self._find_node(profile, 's__Shigella flexneri').coverage)
        self.assertEqual(0.0, self._find_node(profile, 'g__Escherichia').coverage)
        # Total rises only by the residual (2.0 injected - 0.5 reconciled = 1.5)
        total_after = sum([n.coverage for n in profile.breadth_first_iter()])
        self.assertAlmostEqual(total_before + 1.5, total_after)

    # ---- Full reconciliation (--sylph-reconcile) ----

    def test_reconcile_lifts_shared_species_to_sylph(self):
        # E. coli at SingleM 5.0 with 3.0 of genus novel budget; sylph credits
        # 7.0 -> lift to 7.0, drawing the 2.0 increase from genus novel (3 -> 1).
        root = WordNode(None, 'Root')
        root.add_words(['Root'] + self.S_ECOLI.split(';'), 5.0)
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 3.0)
        profile = CondensedCommunityProfile('sample1', root)

        sylph_hits = {_canonical_species_key(self.S_ECOLI): SylphHit(self.S_ECOLI, 7.0)}
        Condenser()._reconcile_with_sylph(profile, sylph_hits, alpha=1.0)

        self.assertEqual(7.0, self._find_node(profile, 's__Escherichia coli').coverage)
        self.assertEqual(1.0, self._find_node(profile, 'g__Escherichia').coverage)

    def test_reconcile_keeps_singlem_when_sylph_is_lower(self):
        # Sylph credits less than SingleM -> SingleM's estimate is kept unchanged
        # (a shared species is never reduced) and novel coverage is untouched.
        root = WordNode(None, 'Root')
        root.add_words(['Root'] + self.S_ECOLI.split(';'), 5.0)
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 3.0)
        profile = CondensedCommunityProfile('sample1', root)

        sylph_hits = {_canonical_species_key(self.S_ECOLI): SylphHit(self.S_ECOLI, 2.0)}
        Condenser()._reconcile_with_sylph(profile, sylph_hits, alpha=1.0)

        self.assertEqual(5.0, self._find_node(profile, 's__Escherichia coli').coverage)
        self.assertEqual(3.0, self._find_node(profile, 'g__Escherichia').coverage)

    def test_reconcile_also_injects_sylph_only_species(self):
        # Reconciliation still adds sylph-only species, like injection.
        root = WordNode(None, 'Root')
        root.add_words(['Root'] + self.S_ECOLI.split(';'), 5.0)
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 3.0)
        profile = CondensedCommunityProfile('sample1', root)

        sylph_hits = {
            _canonical_species_key(self.S_ECOLI): SylphHit(self.S_ECOLI, 5.0),       # equal -> unchanged
            _canonical_species_key(self.S_SHIGELLA): SylphHit(self.S_SHIGELLA, 2.0),  # sylph-only -> injected
        }
        Condenser()._reconcile_with_sylph(profile, sylph_hits, alpha=1.0)

        self.assertEqual(5.0, self._find_node(profile, 's__Escherichia coli').coverage)
        self.assertEqual(2.0, self._find_node(profile, 's__Shigella flexneri').coverage)
        self.assertEqual(1.0, self._find_node(profile, 'g__Escherichia').coverage)

    def test_sylph_profile_read_tsv(self):
        sample_to_hits = SylphProfile.read_tsv(os.path.join(path_to_data, 'small_sylph_profile.tsv'))
        self.assertEqual(['sample1'], list(sample_to_hits.keys()))
        hits = sample_to_hits['sample1']
        self.assertEqual(2, len(hits))
        ecoli_key = _canonical_species_key(self.S_ECOLI)
        self.assertIn(ecoli_key, hits)
        self.assertEqual(9.0, hits[ecoli_key].eff_cov)

    # ---- Joint NNLS deconvolution (--joint) ----

    JOINT_SP1 = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S1'
    JOINT_SP2 = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S2'
    JOINT_GENUS = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G'

    def _joint_otus(self, rows):
        '''rows: list of (coverage, equal_best_list, method). Returns an ArchiveOtuTable.'''
        from singlem.pipe import DIAMOND_ASSIGNMENT_METHOD
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = []
        for i, (coverage, equal_best, method) in enumerate(rows):
            otus.data.append(
                str.split('g{} sample1 seq{}'.format(i, i)) + [1, coverage, '', '', '', '', '', equal_best, method])
        return otus

    def _coverage_of(self, profile, word):
        node = self._find_node(profile, word)
        return node.coverage if node is not None else 0.0

    def test_joint_splits_shared_window_toward_sylph(self):
        from singlem.condense_joint import JointDeconvolver
        # One window shared by S1 and S2; sylph only supports S1.
        otus = self._joint_otus([(5.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD)])
        sylph_hits = {_canonical_species_key(self.JOINT_SP1): SylphHit(self.JOINT_SP1, 5.0)}
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, sylph_hits, alpha=1.0)
        cov = deconv.coverage_by_key
        self.assertGreater(cov[_canonical_species_key(self.JOINT_SP1)], cov.get(_canonical_species_key(self.JOINT_SP2), 0.0))
        self.assertLess(cov.get(_canonical_species_key(self.JOINT_SP2), 0.0), 0.5)

    def test_joint_injects_sylph_only_species(self):
        from singlem.condense_joint import JointDeconvolver
        # No SingleM evidence at all; sylph reports S1 at eff_cov 4. With l1=1,
        # the soft threshold gives a = 4 - l1/2 = 3.5.
        otus = self._joint_otus([])
        sylph_hits = {_canonical_species_key(self.JOINT_SP1): SylphHit(self.JOINT_SP1, 4.0)}
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, sylph_hits, alpha=1.0, l1_penalty=1.0)
        self.assertAlmostEqual(3.5, deconv.coverage_by_key[_canonical_species_key(self.JOINT_SP1)], places=2)

    def test_joint_absent_species_routed_to_novel(self):
        from singlem.condense_joint import JointDeconvolver
        # A window equal-best between species S1 and the genus-level (novel) clade.
        # S1 is absent from sylph, so its coverage is routed to the novel-at-genus leaf.
        otus = self._joint_otus([(5.0, [self.JOINT_SP1, self.JOINT_GENUS], QUERY_BASED_ASSIGNMENT_METHOD)])
        sylph_hits = {}  # S1 not reported by sylph
        deconv = JointDeconvolver()
        profile = deconv.solve('sample1', otus, sylph_hits, alpha=1.0, l1_penalty=0.5, absence_weight=10.0)
        cov = deconv.coverage_by_key
        self.assertLess(cov[_canonical_species_key(self.JOINT_SP1)], 0.5)
        self.assertGreater(cov[_canonical_species_key(self.JOINT_GENUS)], 3.0)
        # In the tree the novel coverage sits on the genus node, not a species leaf.
        self.assertGreater(self._coverage_of(profile, 'g__G'), 3.0)

    def test_joint_min_markers_suppresses_few_marker_species(self):
        from singlem.condense_joint import JointDeconvolver
        # S2 is absent from sylph and uniquely resolved by only one marker (the
        # rest of its support is shared with the sylph-detected S1).
        otus = self._joint_otus([
            (5.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),  # shared marker
            (8.0, [self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),                  # one unique marker
        ])
        sylph_hits = {_canonical_species_key(self.JOINT_SP1): SylphHit(self.JOINT_SP1, 5.0)}
        s2 = _canonical_species_key(self.JOINT_SP2)

        # Default min_markers=3: one unique marker is insufficient -> zeroed.
        d3 = JointDeconvolver()
        d3.solve('s', otus, sylph_hits, alpha=1.0, min_markers=3)
        self.assertEqual(0.0, d3.coverage_by_key.get(s2, 0.0))

        # min_markers=1: the single unique marker now suffices -> retained.
        d1 = JointDeconvolver()
        d1.solve('s', otus, sylph_hits, alpha=1.0, min_markers=1)
        self.assertGreater(d1.coverage_by_key.get(s2, 0.0), 0.0)

    def test_joint_alpha_variable_projection(self):
        from singlem.condense_joint import JointDeconvolver
        # Three species each with a unique window at SingleM coverage 10 and sylph
        # eff_cov 5 -> variable projection should recover alpha ~ 0.5.
        species = [self.JOINT_SP1, self.JOINT_SP2,
                   'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S3']
        otus = self._joint_otus([(10.0, [s], QUERY_BASED_ASSIGNMENT_METHOD) for s in species])
        sylph_hits = {_canonical_species_key(s): SylphHit(s, 5.0) for s in species}
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, sylph_hits, alpha=None, l1_penalty=0.0)
        self.assertAlmostEqual(0.5, deconv.fitted_alpha, places=2)

    # End-to-end Regime 3 test. Reads the mock-metagenome outputs produced by
    # test/data/condense/regime3/Snakefile (run that workflow first) into
    # condense and confirms both the high-coverage genome (recovered by SingleM)
    # and the low-coverage, sylph-only genome (injected by Regime 3) appear in
    # the taxonomic profile.
    REGIME3_METAPACKAGE = '/work/microbiome/db/singlem/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb'

    def test_condense_regime3_mock_metagenome(self):
        import tempfile
        regime3_output = os.path.join(path_to_data, 'regime3', 'output')
        archive = os.path.join(regime3_output, 'archive.json')
        sylph = os.path.join(regime3_output, 'sylph_annotated.tsv')
        for required in (archive, sylph, self.REGIME3_METAPACKAGE):
            if not os.path.exists(required):
                self.skipTest("Regime 3 input not present ({}); run test/data/condense/regime3/Snakefile first".format(required))

        with tempfile.NamedTemporaryFile(suffix='.profile.tsv', mode='w') as profile:
            extern.run("{} condense --input-archive-otu-table {} --metapackage {} "
                "--sylph-profile {} -p {}".format(
                    path_to_script, archive, self.REGIME3_METAPACKAGE, sylph, profile.name))
            with open(profile.name) as f:
                output = f.read()

        # High-coverage genome (10x) detected directly by SingleM.
        self.assertIn('s__Methanobacterium_B sp000744455', output)
        # Low-coverage genome (0.5x): below SingleM's marker sensitivity, so
        # recovered only via sylph and injected by Regime 3.
        self.assertIn('s__Methanobacterium_B lacus', output)

    def test_condense_joint_mock_metagenome(self):
        import tempfile
        regime3_output = os.path.join(path_to_data, 'regime3', 'output')
        archive = os.path.join(regime3_output, 'archive.json')
        sylph = os.path.join(regime3_output, 'sylph_annotated.tsv')
        for required in (archive, sylph, self.REGIME3_METAPACKAGE):
            if not os.path.exists(required):
                self.skipTest("Joint input not present ({}); run test/data/condense/regime3/Snakefile first".format(required))

        with tempfile.NamedTemporaryFile(suffix='.profile.tsv', mode='w') as profile:
            extern.run("{} condense --joint --input-archive-otu-table {} --metapackage {} "
                "--sylph-profile {} -p {}".format(
                    path_to_script, archive, self.REGIME3_METAPACKAGE, sylph, profile.name))
            with open(profile.name) as f:
                output = f.read()

        # Both genomes resolved by the joint deconvolution.
        self.assertIn('s__Methanobacterium_B sp000744455', output)
        self.assertIn('s__Methanobacterium_B lacus', output)

if __name__ == "__main__":
    import logging
    # logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    unittest.main()