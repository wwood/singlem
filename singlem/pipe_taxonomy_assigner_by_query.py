import logging
import os

from .querier import Querier, QueryInputSequence
from .sequence_database import *
from .taxonomy import TaxonomyUtils

class PipeTaxonomyAssignerByQuery:
    def _prepare_query_sequences(self, extracted_reads):
        """ Return at iterable of query sequences to be fed into the querier. """
        spkg_to_queries = {}
        if extracted_reads.analysing_pairs:
            aligned_seqs_to_package_and_sample_name = [{},{}]
        else:
            aligned_seqs_to_package_and_sample_name = [{}]
        for (spkg, extracted_readsets) in extracted_reads.each_package_wise():
            spkg_key = spkg.base_directory()
            if spkg_key not in spkg_to_queries:
                if extracted_reads.analysing_pairs:
                    spkg_to_queries[spkg_key] = [[],[]]
                else:
                    spkg_to_queries[spkg_key] = []
            for maybe_paired_readset in extracted_readsets:
                if extracted_reads.analysing_pairs:
                    for (pair_index, readset) in enumerate(maybe_paired_readset):
                        for read in readset.unknown_sequences:
                            spkg_to_queries[spkg_key][pair_index].append(QueryInputSequence(
                                read.name, read.aligned_sequence, spkg.graftm_package_basename()))
                            # TODO: don't store many copies of the sample name
                            # TODO: Does this fail when there are 2 reads with the same name from different samples?
                            # Always take the sample name as the first of the pair's sample name
                            aligned_seqs_to_package_and_sample_name[pair_index][read.aligned_sequence] = \
                                (spkg, maybe_paired_readset[0].sample_name)
                else:
                    for read in maybe_paired_readset.unknown_sequences:
                        spkg_to_queries[spkg_key].append(QueryInputSequence(read.name, read.aligned_sequence, spkg.graftm_package_basename()))
                        # TODO: don't store many copies of the sample name
                        # TODO: Does this fail when there are 2 reads with the same name from different samples?
                        aligned_seqs_to_package_and_sample_name[0][read.aligned_sequence] = (spkg, maybe_paired_readset.sample_name)

        return spkg_to_queries, aligned_seqs_to_package_and_sample_name

    def assign_taxonomy(self, extracted_reads, assignment_singlem_db, method):
        # query_by_sequence_similarity_with_annoy
        # def query_by_sequence_similarity_with_annoy(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours=None, limit_per_sequence=None):
        
        logging.debug("Preparing query sequences...")
        spkg_queries, aligned_seqs_to_package_and_sample_name = self._prepare_query_sequences(extracted_reads)

        sdb = SequenceDatabase.acquire(assignment_singlem_db)

        logging.debug("Querying...")
        querier = Querier()
        analysing_pairs = extracted_reads.analysing_pairs
        if analysing_pairs:
            final_result = [{},{}]
        else:
            final_result = [{}]

        # TODO: test different values of max_search_nearest_neighbours
        def process_hits_batch(spkg_key, current_hits, pair_index):
            # Get LCA of taxonomy of best hits
            hit_taxonomies = list([h.subject.taxonomy for h in current_hits])
            # We want the final result to be a hash of spkg to sample name to hash of sequence name to taxonomies list
            for hit in current_hits:
                sample_name = aligned_seqs_to_package_and_sample_name[pair_index][hit.query.sequence][1]
                if spkg_key not in final_result[pair_index]:
                    final_result[pair_index][spkg_key] = {}
                if sample_name not in final_result[pair_index][spkg_key]:
                    final_result[pair_index][spkg_key][sample_name] = {}
                final_result[pair_index][spkg_key][sample_name][hit.query.name] = hit_taxonomies

        def query_single_set(queries, pair_index):
            last_query = None
            last_hits = []

            if len(queries) > 0:
                for hit in querier.query_with_queries(queries, sdb, 3, method, SequenceDatabase.NUCLEOTIDE_TYPE, 1, None, False, None):
                    # hit has (query, subject, divergence)
                    # subject has .taxonomy
                    if last_query != hit.query.name:
                        if last_query is not None:
                            # Process last hits batch
                            process_hits_batch(spkg_key, last_hits, pair_index)
                        last_query = hit.query.name
                        last_hits = [hit]
                    else:
                        last_hits.append(hit)
                # Process the last hit
                if last_query is not None:
                    process_hits_batch(spkg_key, last_hits, pair_index)

        for (spkg_key, queries) in spkg_queries.items():
            if analysing_pairs:
                query_single_set(queries[0], 0)
                query_single_set(queries[1], 1)
            else:
                query_single_set(queries, 0)

        if analysing_pairs:
            return QueryTaxonomicAssignmentResult(final_result, analysing_pairs)
        else:
            return QueryTaxonomicAssignmentResult(final_result[0], analysing_pairs)

class QueryTaxonomicAssignmentResult:
    def __init__(self, spkg_to_sample_to_name_to_taxonomies, analysing_pairs):
        self._analysing_pairs = analysing_pairs
        if analysing_pairs:
            # incoming is 2 hashes, we want one hash that ends in a pair of name to taxonomy hashes
            self._spkg_to_sample_to_name_to_taxonomies = {}
            for (spkg_key, sample_to_name_to_taxonomies) in spkg_to_sample_to_name_to_taxonomies[0].items():
                self._spkg_to_sample_to_name_to_taxonomies[spkg_key] = {}
                for (sample_name, name_to_taxonomies) in sample_to_name_to_taxonomies.items():
                    # Add an empty hash for the second readset, in case there are no hits for it.
                    self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name] = [name_to_taxonomies, {}]
            for (spkg_key, sample_to_name_to_taxonomies) in spkg_to_sample_to_name_to_taxonomies[1].items():
                if spkg_key not in self._spkg_to_sample_to_name_to_taxonomies:
                    self._spkg_to_sample_to_name_to_taxonomies[spkg_key] = {}
                for (sample_name, name_to_taxonomies) in sample_to_name_to_taxonomies.items():
                    if sample_name not in self._spkg_to_sample_to_name_to_taxonomies[spkg_key]:
                        # If only read2 has a hit, we need to add an empty read1 hash
                        self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name] = [{}, name_to_taxonomies]
                    else:
                        self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name][1] = name_to_taxonomies
        else:
            self._spkg_to_sample_to_name_to_taxonomies = spkg_to_sample_to_name_to_taxonomies

    def _lca_taxonomy(self, taxonomy_strings):
        lca = TaxonomyUtils.lca_taxonomy_of_strings(taxonomy_strings)
        if lca == []:
            return 'Root'
        else:
            if lca.startswith('Root'):
                return lca
            else:
                return 'Root; '+lca

    def get_best_hits(self, singlem_package, sample_name):
        """ Return dict of read name to LCA of taxonomic hits (or list of 2 dicts for paired reads) """
        spkg_key = singlem_package.base_directory()
        if spkg_key not in self._spkg_to_sample_to_name_to_taxonomies:
            if self._analysing_pairs:
                return [{}, {}]
            else:
                return {}
        if self._analysing_pairs:
            return [{
                name: self._lca_taxonomy(taxonomies)
                for (name, taxonomies) in name_to_taxonomies.items()
            } for name_to_taxonomies in self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name]]
        else:
            return {
                name: self._lca_taxonomy(taxonomies)
                for (name, taxonomies) in self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name].items()
            }

    def get_equal_best_hits(self, singlem_package, sample_name):
        """ Return dict of read name to list of all equal-best taxonomic hits (or list of 2 dicts for paired reads) """
        spkg_key = singlem_package.base_directory()
        if spkg_key not in self._spkg_to_sample_to_name_to_taxonomies:
            return {}
        # Add Root to be consistent with the lca
        if self._analysing_pairs:
            return [{
                name: ['Root; '+tax for tax in taxonomies]
                for (name, taxonomies) in name_to_taxonomies.items()
            } for name_to_taxonomies in self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name]]
        else:
            return {
                name: ['Root; '+tax for tax in taxonomies]
                for (name, taxonomies) in self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name].items()
            }

    def is_assigned_taxonomy(self, singlem_package, sample_name, sequence_name, pair_index):
        spkg_key = singlem_package.base_directory()
        if spkg_key not in self._spkg_to_sample_to_name_to_taxonomies:
            # When no hits are found at all from that pkg
            return False
        if pair_index is None:
            # single-ended case
            return sequence_name in self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name]
        else:
            # paired-ended case
            return sequence_name in self._spkg_to_sample_to_name_to_taxonomies[spkg_key][sample_name][pair_index]