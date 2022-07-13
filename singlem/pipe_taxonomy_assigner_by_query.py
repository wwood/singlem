import logging
import os

from .querier import Querier, QueryInputSequence
from .sequence_database import *

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
            lca = self._lca_taxonomy([h.subject.taxonomy for h in current_hits])
            # We want the final result to be a hash of spkg to sample name to hash of sequence name to taxonomy
            for hit in current_hits:
                sample_name = aligned_seqs_to_package_and_sample_name[pair_index][hit.query.sequence][1]
                if spkg_key not in final_result[pair_index]:
                    final_result[pair_index][spkg_key] = {}
                if sample_name not in final_result[pair_index][spkg_key]:
                    final_result[pair_index][spkg_key][sample_name] = {}
                final_result[pair_index][spkg_key][sample_name][hit.query.name] = lca

        def query_single_set(queries, pair_index):
            last_query = None
            last_hits = []

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

    def _lca_taxonomy(self, taxonomy_strings):
        hit_taxonomies = list([list([ta.strip() for ta in t.split(';')]) for t in taxonomy_strings])
        lca = hit_taxonomies[0]
        for taxonomy in hit_taxonomies[1:]:
            if len(taxonomy) < len(lca):
                lca = lca[:len(taxonomy)]
            for i, tax in enumerate(taxonomy):
                if i >= len(lca) or tax != lca[i]:
                    lca = lca[:i]
                    break
        if lca == []:
            return 'Root'
        else:
            return 'Root; '+'; '.join(lca)

class QueryTaxonomicAssignmentResult:
    def __init__(self, spkg_to_sample_to_name_to_taxonomy, analysing_pairs):
        self._analysing_pairs = analysing_pairs
        if analysing_pairs:
            # incoming is 2 hashes, we want one hash that ends in a pair of name to taxonomy hashes
            self._spkg_to_sample_to_name_to_taxonomy = {}
            for (spkg_key, sample_to_name_to_taxonomy) in spkg_to_sample_to_name_to_taxonomy[0].items():
                self._spkg_to_sample_to_name_to_taxonomy[spkg_key] = {}
                for (sample_name, name_to_taxonomy) in sample_to_name_to_taxonomy.items():
                    self._spkg_to_sample_to_name_to_taxonomy[spkg_key][sample_name] = [name_to_taxonomy]
            for (spkg_key, sample_to_name_to_taxonomy) in spkg_to_sample_to_name_to_taxonomy[1].items():
                for (sample_name, name_to_taxonomy) in sample_to_name_to_taxonomy.items():
                    self._spkg_to_sample_to_name_to_taxonomy[spkg_key][sample_name].append(name_to_taxonomy)
        else:
            self._spkg_to_sample_to_name_to_taxonomy = spkg_to_sample_to_name_to_taxonomy

    def get_best_hits(self, singlem_package, sample_name):
        # if self._analysing_pairs:
        #     return self._spkg_to_sample_to_name_to_taxonomy[singlem_package.base_directory()][sample_name]
        # else:
        spkg_key = singlem_package.base_directory()
        if spkg_key not in self._spkg_to_sample_to_name_to_taxonomy:
            return {}
        return self._spkg_to_sample_to_name_to_taxonomy[singlem_package.base_directory()][sample_name]

    def is_assigned_taxonomy(self, singlem_package, sample_name, sequence_name, pair_index):
        spkg_key = singlem_package.base_directory()
        if spkg_key not in self._spkg_to_sample_to_name_to_taxonomy:
            # When no hits are found at all from that pkg
            return False
        if pair_index is None:
            # single-ended case
            return sequence_name in self._spkg_to_sample_to_name_to_taxonomy[spkg_key][sample_name]
        else:
            # paired-ended case
            return sequence_name in self._spkg_to_sample_to_name_to_taxonomy[spkg_key][sample_name][pair_index]