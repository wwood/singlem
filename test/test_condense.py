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
from singlem.condense  import Condenser, WordNode, CondensedCommunityProfile, WeebillHit, WeebillProfile
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

    # ---- Regime 3: weebill-only species injection ----

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
        weebill = {'a': WeebillHit('a', 5.0), 'b': WeebillHit('b', 10.0), 'c': WeebillHit('c', 15.0)}
        self.assertAlmostEqual(2.0, Condenser()._fit_alpha(singlem, weebill))

    def test_fit_alpha_too_few_anchors(self):
        # Only two species at >= 10x coverage => default to 1.0
        singlem = {'a': 10.0, 'b': 20.0, 'c': 3.0}
        weebill = {'a': WeebillHit('a', 5.0), 'b': WeebillHit('b', 10.0), 'c': WeebillHit('c', 1.5)}
        self.assertEqual(1.0, Condenser()._fit_alpha(singlem, weebill))

    def _find_node(self, profile, word):
        for node in profile.breadth_first_iter():
            if node.word == word:
                return node
        return None

    def test_inject_weebill_only_species_reconciles_with_genus(self):
        root = WordNode(None, 'Root')
        root.add_words(['Root'] + self.S_ECOLI.split(';'), 5.0)
        # Genus-level novel coverage that injection should draw down
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 3.0)
        profile = CondensedCommunityProfile('sample1', root)

        weebill_hits = {
            _canonical_species_key(self.S_ECOLI): WeebillHit(self.S_ECOLI, 9.0),  # already present
            _canonical_species_key(self.S_SHIGELLA): WeebillHit(self.S_SHIGELLA, 2.0),  # weebill-only
        }
        Condenser()._inject_weebill_only_species(profile, weebill_hits, alpha=1.0)

        # E. coli already in profile, untouched
        self.assertEqual(5.0, self._find_node(profile, 's__Escherichia coli').coverage)
        # Shigella injected at alpha*eff_cov = 2.0
        self.assertEqual(2.0, self._find_node(profile, 's__Shigella flexneri').coverage)
        # Genus novel coverage drawn down from 3.0 to 1.0 (residual 0, no new total)
        self.assertEqual(1.0, self._find_node(profile, 'g__Escherichia').coverage)

    def test_inject_weebill_only_species_residual_is_new_coverage(self):
        root = WordNode(None, 'Root')
        # Only 0.5 of genus novel budget available
        root.add_words(['Root'] + self.G_ECOLI.split(';'), 0.5)
        profile = CondensedCommunityProfile('sample1', root)
        total_before = sum([n.coverage for n in profile.breadth_first_iter()])

        weebill_hits = {_canonical_species_key(self.S_SHIGELLA): WeebillHit(self.S_SHIGELLA, 2.0)}
        Condenser()._inject_weebill_only_species(profile, weebill_hits, alpha=1.0)

        self.assertEqual(2.0, self._find_node(profile, 's__Shigella flexneri').coverage)
        self.assertEqual(0.0, self._find_node(profile, 'g__Escherichia').coverage)
        # Total rises only by the residual (2.0 injected - 0.5 reconciled = 1.5)
        total_after = sum([n.coverage for n in profile.breadth_first_iter()])
        self.assertAlmostEqual(total_before + 1.5, total_after)

    def test_weebill_profile_read_tsv(self):
        sample_to_hits = WeebillProfile.read_tsv(os.path.join(path_to_data, 'small_weebill_profile.tsv'))
        self.assertEqual(['sample1'], list(sample_to_hits.keys()))
        hits = sample_to_hits['sample1']
        self.assertEqual(2, len(hits))
        ecoli_key = _canonical_species_key(self.S_ECOLI)
        self.assertIn(ecoli_key, hits)
        self.assertEqual(9.0, hits[ecoli_key].eff_cov)

    def test_weebill_profile_read_tsv_true_cov(self):
        # weebill -u (--estimate-unknown) reports True_cov in place of Eff_cov: a
        # coverage corrected for the unknown sequence fraction, and so already on
        # SingleM's scale. The column name is the only record of which was run, so
        # it is preserved through annotation and both are accepted here.
        sample_to_hits = WeebillProfile.read_tsv(
            os.path.join(path_to_data, 'small_weebill_profile_true_cov.tsv'))
        hits = sample_to_hits['sample1']
        self.assertEqual(2, len(hits))
        self.assertEqual(9.0, hits[_canonical_species_key(self.S_ECOLI)].eff_cov)

    def test_weebill_profile_is_unknown_corrected(self):
        # The column name is the only record of whether weebill was run with -u, and
        # condense keys alpha off it: True_cov is already on SingleM's coverage
        # scale, Eff_cov is a fixed fraction of it and must be calibrated.
        self.assertTrue(WeebillProfile.is_unknown_corrected(
            os.path.join(path_to_data, 'small_weebill_profile_true_cov.tsv')))
        self.assertFalse(WeebillProfile.is_unknown_corrected(
            os.path.join(path_to_data, 'small_weebill_profile.tsv')))

    # ---- Joint NNLS deconvolution (--joint) ----

    JOINT_SP1 = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S1'
    JOINT_SP2 = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S2'
    JOINT_SP3 = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S3'
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

    def test_joint_splits_shared_window_toward_weebill(self):
        from singlem.condense_joint import JointDeconvolver
        # One window shared by S1 and S2; weebill only supports S1.
        otus = self._joint_otus([(5.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD)])
        weebill_hits = {_canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 5.0)}
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, weebill_hits, alpha=1.0)
        cov = deconv.coverage_by_key
        self.assertGreater(cov[_canonical_species_key(self.JOINT_SP1)], cov.get(_canonical_species_key(self.JOINT_SP2), 0.0))
        self.assertLess(cov.get(_canonical_species_key(self.JOINT_SP2), 0.0), 0.5)

    def test_joint_injects_weebill_only_species(self):
        from singlem.condense_joint import JointDeconvolver
        # No SingleM evidence at all; weebill reports S1 at eff_cov 4. The soft threshold
        # gives a = eff_cov - l1_penalty / (2 * weebill_weight): the more the model defers
        # to weebill, the less the l1 penalty shrinks a species weebill alone supports.
        otus = self._joint_otus([])
        weebill_hits = {_canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 4.0)}
        for weebill_weight, expected in [(1.0, 3.5), (50.0, 3.99)]:
            deconv = JointDeconvolver()
            deconv.solve('sample1', otus, weebill_hits, alpha=1.0, l1_penalty=1.0,
                         weebill_weight=weebill_weight)
            self.assertAlmostEqual(
                expected, deconv.coverage_by_key[_canonical_species_key(self.JOINT_SP1)], places=2)

    def test_joint_absent_species_routed_to_novel(self):
        from singlem.condense_joint import JointDeconvolver
        # A window equal-best between species S1 and the genus-level (novel) clade.
        # S1 is absent from weebill, so its coverage is routed to the novel-at-genus leaf.
        # The marker floor is disabled because this sample is a single window that is
        # shared between S1 and G, so it is unique evidence for neither and the floor
        # (which applies to novel columns as well as species) would zero both. The floor
        # is tested separately below; what is under test here is where the coverage is
        # routed, not how much evidence there has to be.
        otus = self._joint_otus([(5.0, [self.JOINT_SP1, self.JOINT_GENUS], QUERY_BASED_ASSIGNMENT_METHOD)])
        weebill_hits = {}  # S1 not reported by weebill
        deconv = JointDeconvolver()
        profile = deconv.solve('sample1', otus, weebill_hits, alpha=1.0, l1_penalty=0.5, absence_weight=10.0,
                               min_markers=0)
        cov = deconv.coverage_by_key
        self.assertLess(cov[_canonical_species_key(self.JOINT_SP1)], 0.5)
        self.assertGreater(cov[_canonical_species_key(self.JOINT_GENUS)], 3.0)
        # In the tree the novel coverage sits on the genus node, not a species leaf.
        self.assertGreater(self._coverage_of(profile, 'g__G'), 3.0)

    def test_joint_min_markers_suppresses_few_marker_species(self):
        from singlem.condense_joint import JointDeconvolver
        # S2 is absent from weebill and uniquely resolved by only one marker (the
        # rest of its support is shared with the weebill-detected S1).
        otus = self._joint_otus([
            (5.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),  # shared marker
            (8.0, [self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),                  # one unique marker
        ])
        weebill_hits = {_canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 5.0)}
        s2 = _canonical_species_key(self.JOINT_SP2)

        # Default min_markers=3: one unique marker is insufficient -> zeroed.
        d3 = JointDeconvolver()
        d3.solve('s', otus, weebill_hits, alpha=1.0, min_markers=3)
        self.assertEqual(0.0, d3.coverage_by_key.get(s2, 0.0))

        # min_markers=1: the single unique marker now suffices -> retained.
        d1 = JointDeconvolver()
        d1.solve('s', otus, weebill_hits, alpha=1.0, min_markers=1,
                 min_singlem_coverage=0.0)
        self.assertGreater(d1.coverage_by_key.get(s2, 0.0), 0.0)

    # Genus-level novel columns in three different phyla; their only shared rank is the
    # domain, so a rollup settles them there.
    ROLLUP_G1 = 'd__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1'
    ROLLUP_G2 = 'd__Bacteria;p__P2;c__C2;o__O2;f__F2;g__G2'
    ROLLUP_G3 = 'd__Bacteria;p__P3;c__C3;o__O3;f__F3;g__G3'

    def test_joint_scattered_novel_markers_roll_up_to_common_ancestor(self):
        from singlem.condense_joint import JointDeconvolver
        # A phylogenetically novel organism with no database representative matches a
        # different database relative on every marker, so its signal scatters into several
        # genus-level novel columns of one marker each. None clears the min_markers floor
        # alone, so without rollup all are zeroed and the organism vanishes. Rolled up,
        # their shared evidence of novelty settles on their common ancestor -- here, since
        # the three genera lie in different phyla, the domain -- and is reported there. The
        # organism is seen on all three of its domain's markers, so its zero-padded trimmed
        # mean is its full ~10x.
        otus = self._joint_otus([
            (10.0, [self.ROLLUP_G1], QUERY_BASED_ASSIGNMENT_METHOD),
            (10.0, [self.ROLLUP_G2], QUERY_BASED_ASSIGNMENT_METHOD),
            (10.0, [self.ROLLUP_G3], QUERY_BASED_ASSIGNMENT_METHOD),
        ])
        deconv = JointDeconvolver()
        profile = deconv.solve('sample1', otus, {}, domain_marker_counts={'Bacteria': 3},
                               alpha=1.0, min_markers=3)
        cov = deconv.coverage_by_key
        for g in (self.ROLLUP_G1, self.ROLLUP_G2, self.ROLLUP_G3):
            self.assertEqual(0.0, cov.get(_canonical_species_key(g), 0.0))
        self.assertGreater(cov.get('d__Bacteria', 0.0), 8.0)
        # The novelty sits on the domain node of the tree, not any genus leaf.
        self.assertGreater(self._coverage_of(profile, 'd__Bacteria'), 8.0)

    def test_joint_rollup_misassignment_sink_is_held_near_zero(self):
        from singlem.condense_joint import JointDeconvolver
        # A clade with nothing truly novel in it still collects scattered genus columns --
        # an abundant organism's reads mis-hitting foreign markers and tying across the
        # clade -- and once rolled up the clade is their sole bidder. But a misassignment
        # sink is seen on only a handful of its domain's markers, even when those carry
        # enormous coverage (here one 1000x spike). Ungoverned the ancestor would fit the
        # mean of its detected markers and fabricate novelty; held instead to the same
        # zero-padded trimmed mean the standard condense uses -- three markers out of a
        # thirty-marker domain, the spike trimmed -- it stays near zero.
        otus = self._joint_otus([
            (1000.0, [self.ROLLUP_G1], QUERY_BASED_ASSIGNMENT_METHOD),  # misassignment spike
            (2.0, [self.ROLLUP_G2], QUERY_BASED_ASSIGNMENT_METHOD),
            (2.0, [self.ROLLUP_G3], QUERY_BASED_ASSIGNMENT_METHOD),
        ])
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, {}, domain_marker_counts={'Bacteria': 30},
                     alpha=1.0, min_markers=3)
        self.assertLess(deconv.coverage_by_key.get('d__Bacteria', 0.0), 1.0)

    def test_joint_species_level_tie_is_evidence_for_the_novel_clade(self):
        from singlem.condense_joint import JointDeconvolver
        # A novel organism in a genus whose species are all in the database need never
        # produce a window that SingleM resolves to the genus: its window can tie across
        # the DB species on every single marker. The tie resolves the read no deeper than
        # the genus, so it is evidence for the genus's novel column, and the marker floor
        # must accept it as such. Otherwise the novel column earns zero markers, is fixed
        # to zero before the optimiser runs, and the organism's coverage disappears from
        # the profile rather than being reported as novel.
        rows = []
        for _ in range(3):
            rows.append((30.0, [self.JOINT_SP1], QUERY_BASED_ASSIGNMENT_METHOD))
            rows.append((100.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD))
        # _joint_otus gives each row its own marker, so pair them up onto three markers.
        otus = self._joint_otus(rows)
        for i, otu in enumerate(otus.data):
            otu[0] = 'g{}'.format(i // 2)
        weebill_hits = {_canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 30.0)}

        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, weebill_hits, domain_marker_counts={'Bacteria': 3},
                     alpha=1.0, min_markers=3)
        cov = deconv.coverage_by_key
        # The novel organism is recovered at close to its true 100x, and S1 is not
        # inflated by it: weebill pins S1 at 30x and the tie's coverage goes to the genus.
        self.assertGreater(cov[_canonical_species_key(self.JOINT_GENUS)], 90.0)
        self.assertAlmostEqual(30.0, cov[_canonical_species_key(self.JOINT_SP1)], delta=1.0)
        # S2, which weebill does not report and which no window resolves to alone, stays out.
        self.assertEqual(0.0, cov.get(_canonical_species_key(self.JOINT_SP2), 0.0))

    def test_joint_diverged_window_of_known_species_is_not_novelty(self):
        from singlem.condense_joint import JointDeconvolver
        # S1 is abundant and confirmed by weebill, and has its own window on m0 and m1.
        # On m2 the sample's strain has diverged from the database representative at that
        # one window, so the read ties across congeners *without naming S1 at all*. Those
        # are S1's reads: single-copy markers are universal, so S1 has a window on m2, and
        # if it is not in a row of its own it is in this one. Unless S1 may claim them,
        # the genus's novel column is the only bidder and ordinary strain variation in a
        # known species is reported as a novel congener.
        otus = self._joint_otus([
            (10.0, [self.JOINT_SP1], QUERY_BASED_ASSIGNMENT_METHOD),
            (10.0, [self.JOINT_SP1], QUERY_BASED_ASSIGNMENT_METHOD),
            (10.0, [self.JOINT_SP2, self.JOINT_SP3], QUERY_BASED_ASSIGNMENT_METHOD),
        ])
        for i, otu in enumerate(otus.data):
            otu[0] = 'g{}'.format(i)
        weebill_hits = {_canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 10.0)}
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, weebill_hits, domain_marker_counts={'Bacteria': 3},
                     alpha=1.0, min_markers=2)
        cov = deconv.coverage_by_key
        # S1 keeps its coverage and the tie does not become a novel organism.
        self.assertAlmostEqual(10.0, cov[_canonical_species_key(self.JOINT_SP1)], delta=1.0)
        self.assertLess(cov.get(_canonical_species_key(self.JOINT_GENUS), 0.0), 3.0)

    def test_joint_weebill_ceiling_stops_rare_species_absorbing_clade_coverage(self):
        from singlem.condense_joint import JointDeconvolver
        # S2 is reported by weebill at 0.2x -- it is barely there -- but it is a candidate
        # in a genus carrying a great deal of ambiguous coverage. The weebill rows alone are
        # a quadratic on an absolute residual, so enough high-coverage SingleM rows can
        # drag S2 far above what weebill says; the ceiling is what stops the rarest species
        # in a clade from becoming the cheapest home for the clade's ambiguous coverage.
        otus = self._joint_otus([
            (100.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),
            (100.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),
            (100.0, [self.JOINT_SP1, self.JOINT_SP2], QUERY_BASED_ASSIGNMENT_METHOD),
        ])
        for i, otu in enumerate(otus.data):
            otu[0] = 'g{}'.format(i)
        weebill_hits = {
            _canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 100.0),
            _canonical_species_key(self.JOINT_SP2): WeebillHit(self.JOINT_SP2, 0.2),
        }
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, weebill_hits, domain_marker_counts={'Bacteria': 3},
                     alpha=1.0, min_markers=0)
        cov = deconv.coverage_by_key
        # S2 stays near weebill's 0.2x (the ceiling allows slack, not a free hand), and the
        # ambiguous coverage goes to S1, which weebill says is the abundant one.
        self.assertLess(cov[_canonical_species_key(self.JOINT_SP2)], 2.0)
        self.assertGreater(cov[_canonical_species_key(self.JOINT_SP1)], 90.0)

    def test_joint_filters_low_coverage_singlem_only_genus(self):
        from singlem.condense_joint import JointDeconvolver
        otus = self._joint_otus([
            (0.2, [self.JOINT_GENUS], QUERY_BASED_ASSIGNMENT_METHOD),
        ])
        deconv = JointDeconvolver()
        profile = deconv.solve(
            'sample1', otus, {}, alpha=1.0, l1_penalty=0.0,
            min_singlem_coverage=0.35)
        self.assertEqual(0.0, deconv.coverage_by_key[_canonical_species_key(self.JOINT_GENUS)])
        self.assertEqual(0.0, self._coverage_of(profile, 'g__G'))

    def test_joint_retains_low_coverage_weebill_species(self):
        from singlem.condense_joint import JointDeconvolver
        otus = self._joint_otus([])
        weebill_hits = {_canonical_species_key(self.JOINT_SP1): WeebillHit(self.JOINT_SP1, 0.2)}
        deconv = JointDeconvolver()
        profile = deconv.solve(
            'sample1', otus, weebill_hits, alpha=1.0, l1_penalty=0.0,
            min_singlem_coverage=0.35)
        self.assertAlmostEqual(0.2, deconv.coverage_by_key[_canonical_species_key(self.JOINT_SP1)], places=3)
        self.assertAlmostEqual(0.2, self._coverage_of(profile, 's__S1'), places=3)

    def test_joint_alpha_variable_projection(self):
        from singlem.condense_joint import JointDeconvolver
        # Three species each with a unique window at SingleM coverage 10 and weebill
        # eff_cov 5 -> variable projection should recover alpha ~ 0.5.
        species = [self.JOINT_SP1, self.JOINT_SP2,
                   'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S3']
        otus = self._joint_otus([(10.0, [s], QUERY_BASED_ASSIGNMENT_METHOD) for s in species])
        weebill_hits = {_canonical_species_key(s): WeebillHit(s, 5.0) for s in species}
        deconv = JointDeconvolver()
        deconv.solve('sample1', otus, weebill_hits, alpha=None, l1_penalty=0.0)
        self.assertAlmostEqual(0.5, deconv.fitted_alpha, places=2)

    def test_joint_weebill_leverage_does_not_depend_on_alpha(self):
        from singlem.condense_joint import JointDeconvolver
        # S1's markers all read 10x but weebill puts it at 2x in SingleM's units. How far
        # the fit is pulled off the markers toward weebill is set by weebill_weight, and must
        # not also depend on alpha -- the residual is taken as (e/alpha - a), so a
        # four-fold change in alpha with e scaled to match leaves the problem unchanged.
        # Taking it as (e - alpha*a) instead would scale the weebill term by alpha^2 and
        # quarter weebill's influence here, which is the wrong direction: alpha is smallest
        # in the low-coverage samples where the markers deserve the least trust.
        otus = self._joint_otus(
            [(10.0, [self.JOINT_SP1], QUERY_BASED_ASSIGNMENT_METHOD) for _ in range(10)])
        s1 = _canonical_species_key(self.JOINT_SP1)
        fitted = []
        for alpha in (1.0, 0.25):
            deconv = JointDeconvolver()
            deconv.solve('sample1', otus, {s1: WeebillHit(self.JOINT_SP1, 2.0 * alpha)},
                         alpha=alpha, domain_marker_counts={'Bacteria': 10})
            fitted.append(deconv.coverage_by_key[s1])
        self.assertAlmostEqual(fitted[0], fitted[1], places=3)
        # And weebill, weighted 150 against 10 markers, is what the fit mostly follows.
        self.assertLess(fitted[0], 4.0)

    # End-to-end Regime 3 test. Reads the mock-metagenome outputs produced by
    # test/data/condense/regime3/Snakefile (run that workflow first) into
    # condense and confirms both the high-coverage genome (recovered by SingleM)
    # and the low-coverage, weebill-only genome (injected by Regime 3) appear in
    # the taxonomic profile.
    REGIME3_METAPACKAGE = '/work/microbiome/db/singlem/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb'

    def test_condense_regime3_mock_metagenome(self):
        import tempfile
        regime3_output = os.path.join(path_to_data, 'regime3', 'output')
        archive = os.path.join(regime3_output, 'archive.json')
        weebill = os.path.join(regime3_output, 'weebill_annotated.tsv')
        for required in (archive, weebill, self.REGIME3_METAPACKAGE):
            if not os.path.exists(required):
                self.skipTest("Regime 3 input not present ({}); run test/data/condense/regime3/Snakefile first".format(required))

        with tempfile.NamedTemporaryFile(suffix='.profile.tsv', mode='w') as profile:
            extern.run("{} condense --input-archive-otu-table {} --metapackage {} "
                "--weebill-profile {} -p {}".format(
                    path_to_script, archive, self.REGIME3_METAPACKAGE, weebill, profile.name))
            with open(profile.name) as f:
                output = f.read()

        # High-coverage genome (10x) detected directly by SingleM.
        self.assertIn('s__Methanobacterium_B sp000744455', output)
        # Low-coverage genome (0.5x): below SingleM's marker sensitivity, so
        # recovered only via weebill and injected by Regime 3.
        self.assertIn('s__Methanobacterium_B lacus', output)

    def test_condense_joint_mock_metagenome(self):
        import tempfile
        regime3_output = os.path.join(path_to_data, 'regime3', 'output')
        archive = os.path.join(regime3_output, 'archive.json')
        weebill = os.path.join(regime3_output, 'weebill_annotated.tsv')
        for required in (archive, weebill, self.REGIME3_METAPACKAGE):
            if not os.path.exists(required):
                self.skipTest("Joint input not present ({}); run test/data/condense/regime3/Snakefile first".format(required))

        with tempfile.NamedTemporaryFile(suffix='.profile.tsv', mode='w') as profile:
            extern.run("{} condense --joint --input-archive-otu-table {} --metapackage {} "
                "--weebill-profile {} -p {}".format(
                    path_to_script, archive, self.REGIME3_METAPACKAGE, weebill, profile.name))
            with open(profile.name) as f:
                output = f.read()

        # Both genomes resolved by the joint deconvolution.
        self.assertIn('s__Methanobacterium_B sp000744455', output)
        self.assertIn('s__Methanobacterium_B lacus', output)

if __name__ == "__main__":
    import logging
    # logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    unittest.main()
