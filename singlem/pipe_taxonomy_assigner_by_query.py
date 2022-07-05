import logging
import os

from .querier import Querier, QueryInputSequence
from .sequence_database import SequenceDatabase

class PipeTaxonomyAssignerByQuery:
    def prepare_query_sequences(self, extracted_reads):
        """ Return at iterable of query sequences to be fed into the querier. """
        spkg_to_queries = {}
        aligned_seqs_to_package_and_sample_name = {}
        for (spkg, extracted_readsets) in extracted_reads.each_package_wise():
            for (readset_i, readset) in enumerate(extracted_readsets):
                for read in readset.unknown_sequences:
                    spkg_key = spkg.base_directory()
                    if spkg_key not in spkg_to_queries:
                        spkg_to_queries[spkg_key] = []
                    spkg_to_queries[spkg_key].append(QueryInputSequence(read.name, read.aligned_sequence, spkg.graftm_package_basename()))
                    # TODO: don't store many copies of the sample name
                    # TODO: Does this fail when there are 2 reads with the same name from different samples?
                    aligned_seqs_to_package_and_sample_name[read.aligned_sequence] = (spkg, readset.sample_name)

        return spkg_to_queries, aligned_seqs_to_package_and_sample_name

    def assign_taxonomy_with_annoy(self, extracted_reads, assignment_singlem_db):
        # query_by_sequence_similarity_with_annoy
        # def query_by_sequence_similarity_with_annoy(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours=None, limit_per_sequence=None):
        
        logging.debug("Preparing query sequences...")
        spkg_queries, aligned_seqs_to_package_and_sample_name = self.prepare_query_sequences(extracted_reads)

        sdb = SequenceDatabase.acquire(assignment_singlem_db)

        logging.debug("Querying...")
        querier = Querier()
        # TODO: test different values of max_search_nearest_neighbours
        last_query = None
        last_hits = []
        final_result = {}
        def process_hits_batch(spkg_key, final_result, current_hits):
            # Get LCA of taxonomy of best hits
            lca = self._lca_taxonomy([h.subject.taxonomy for h in last_hits])
            # We want the final result to be a hash of spkg to sample name to hash of sequence name to taxonomy
            for hit in current_hits:
                sample_name = aligned_seqs_to_package_and_sample_name[hit.query.sequence][1]
                if spkg_key not in final_result:
                    final_result[spkg_key] = {}
                if sample_name not in final_result[spkg_key]:
                    final_result[spkg_key][sample_name] = {}
                final_result[spkg_key][sample_name][hit.query.name] = lca

        for (spkg_key, queries) in spkg_queries.items():
            # for hit in querier.query_by_sequence_similarity_with_annoy(queries, sdb, 3, SequenceDatabase.NUCLEOTIDE_TYPE, 1):
            for hit in querier.query_with_queries(queries, sdb, 3, 'annoy', SequenceDatabase.NUCLEOTIDE_TYPE, 1, None, False, None):
                # hit has (query, subject, divergence)
                # subject has .taxonomy
                if last_query != hit.query.name:
                    if last_query is not None:
                        # Process last hits batch
                        process_hits_batch(spkg_key, final_result, last_hits)
                    last_query = hit.query.name
                    last_hits = [hit]
                else:
                    last_hits.append(hit)
            # Process the last hit
            if last_query is not None:
                process_hits_batch(spkg_key, final_result, last_hits)

        return QueryTaxonomicAssignmentResult(final_result)

    def _lca_taxonomy(self, taxonomy_strings):
        lca = []
        hit_taxonomies = list([list(t.split('; ')) for t in taxonomy_strings])

        for (i, taxon) in enumerate(hit_taxonomies[0]):
            if all([len(h) > i and h[i]==taxon for h in hit_taxonomies]):
                lca.append(taxon)
            else:
                break
        if lca == []:
            return 'Root'
        else:
            return 'Root; '+'; '.join(lca)

class QueryTaxonomicAssignmentResult:
    def __init__(self, spkg_to_sample_to_name_to_taxonomy):
        self._spkg_to_sample_to_name_to_taxonomy = spkg_to_sample_to_name_to_taxonomy

    def get_best_hits(self, singlem_package, sample_name):
        return self._spkg_to_sample_to_name_to_taxonomy[singlem_package.base_directory()][sample_name]