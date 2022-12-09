import tempfile
import logging
import subprocess
import sys
from Bio import pairwise2
import numpy as np
import csv
import pandas as pd
import itertools
import math
from bird_tool_utils import iterable_chunks
from sqlalchemy import select
from sqlalchemy.orm import Session

from .sequence_database import SequenceDatabase
from . import sequence_database
from .singlem_database_models import *
from .sequence_classes import SeqReader
from .query_formatters import SparseResultFormatter
from .otu_table_collection import OtuTableCollection
from .otu_table_entry import OtuTableEntry
from .otu_table import OtuTable
from .otu_table_collection import StreamingOtuTableCollection

class Querier:
    def __init__(self):
        self._query_result_from_db_builder_nucleotide = None
        self._query_result_from_db_builder_protein = None

    def query(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'), min_version=5)
        max_divergence = kwargs.pop('max_divergence')
        output_style = kwargs.pop('output_style')
        query_otu_table = kwargs.pop('query_otu_table')
        num_threads = kwargs.pop('num_threads')
        search_method = kwargs.pop('search_method')
        sequence_type = kwargs.pop('sequence_type')
        # stream_output = kwargs.pop('stream_output')
        max_nearest_neighbours = kwargs.pop('max_nearest_neighbours')
        max_search_nearest_neighbours = kwargs.pop('max_search_nearest_neighbours')
        preload_db = kwargs.pop('preload_db')
        limit_per_sequence = kwargs.pop('limit_per_sequence')

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        queries = self.prepare_query_sequences(query_otu_table, num_threads)

        query_results = self.query_with_queries(
            queries, db, max_divergence, search_method, sequence_type, 
            max_nearest_neighbours=max_nearest_neighbours, max_search_nearest_neighbours=max_search_nearest_neighbours, \
            preload_db=preload_db, limit_per_sequence=limit_per_sequence)
        # Only reason not to stream would be so that the queries are returned in the same order as passed in. eh for now.
        do_stream = True
        if not do_stream:
            query_results = list(query_results)
            logging.info("Printing %i hits" % len(query_results))

        if output_style == 'sparse':
            SparseResultFormatter().write(query_results, sys.stdout, sequence_type, streaming=do_stream)
        elif output_style == None:
            return query_results
        else:
            raise Exception("Programming error")

    def preload_nucleotide_db(self, sdb, marker_id):
        with Session(sdb.sqlalchemy_connection) as conn:
            marker_name = conn.execute(select(Marker.marker).where(Marker.id==marker_id)).fetchone()[0]
            logging.info("Caching nucleotide data for marker {}..".format(marker_name))

            query = select([
                Otu.marker_wise_sequence_id,
                Otu.sequence,
                Otu.sample_name,
                Otu.num_hits,
                Otu.coverage,
                Taxonomy.taxonomy]) \
                    .where(Otu.taxonomy_id == Taxonomy.id) \
                    .where(Otu.marker_id == marker_id)
            result = conn.execute(query)
            d1 = pd.DataFrame(result.fetchall(), 
                columns = ('nucleotides_marker_wise_id','nucleotide_sequence', \
                    'sample_name', 'num_hits', 'coverage', 'taxonomy'))
        return d1

    def preload_protein_db(self, sdb, marker_id):
        with Session(sdb.sqlalchemy_connection) as conn:
            marker_name = conn.execute(select(Marker.marker).where(Marker.id==marker_id)).fetchone()[0]
            logging.info("Caching protein data for marker {}..".format(marker_name))

            query = select([
                ProteinSequence.marker_wise_id,
                Otu.sequence,
                ProteinSequence.protein_sequence,
                Otu.sample_name,
                Otu.num_hits,
                Otu.coverage,
                Taxonomy.taxonomy]) \
                    .where(Otu.taxonomy_id == Taxonomy.id) \
                    .where(Otu.marker_id == marker_id) \
                    .where(NucleotidesProteins.nucleotide_id == Otu.sequence_id) \
                    .where(NucleotidesProteins.protein_id == ProteinSequence.id)
            result = conn.execute(query)
            d1 = pd.DataFrame(result.fetchall(), 
                columns = ('proteins_marker_wise_id','nucleotide_sequence','protein_sequence', \
                'sample_name', 'num_hits', 'coverage', 'taxonomy'))
        return d1




    def prepare_query_sequences(self, otus, num_threads):
        '''return an iterable of QueryInputSequence objects sorted by marker.'''

        # To reduce RAM requirements, sort the input by marker, so that the DBs
        # can be dropped when that marker is finished being queried.
        logging.info("Sorting query sequences by marker ID to save RAM ..")
        sorted_io = tempfile.NamedTemporaryFile(prefix='singlem-query-sorted')
        sorted_path = sorted_io.name
        proc = subprocess.Popen(['bash','-c','sort --parallel={} --buffer-size=20% > {}'.format(num_threads, sorted_path)],
            stdin=subprocess.PIPE,
            stdout=None,
            stderr=subprocess.PIPE,
            universal_newlines=True)
        total_otu_count = 0
        for e in otus:
            total_otu_count += 1
            print("\t".join([
                e.marker,
                e.sequence,
                e.sample_name
            ]), file=proc.stdin)
        proc.stdin.close()
        proc.wait()
        if proc.returncode != 0:
            raise Exception("Sort command returned non-zero exit status %i.\n"\
                "STDERR was: %s" % (
                    proc.returncode, proc.stderr.read()))
        logging.info("Sorted {} OTU queries ..".format(total_otu_count))
        return MarkerSortedQueryInput(sorted_io)

    def query_with_queries(self, queries, sdb, max_divergence, search_method, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours, preload_db, limit_per_sequence):
        if max_divergence == 0 and sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
            if limit_per_sequence != None:
                raise Exception("limit-per-sequence has not been implemented for nucleotide queries with max-divergence 0 yet")
            return self.query_by_sqlite(queries, sdb)
        elif search_method == 'scann-naive':
            return self.query_by_sequence_similarity_with_scann(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, naive=True, preload_db=preload_db, limit_per_sequence=limit_per_sequence)
        elif search_method == 'annoy':
            if preload_db:
                raise NotImplementedError("Preloading the database is not supported for the annoy search method")
            return self.query_by_sequence_similarity_with_annoy(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours=max_search_nearest_neighbours, limit_per_sequence=limit_per_sequence)
        elif search_method == 'nmslib':
            if preload_db:
                raise NotImplementedError("Preloading the database is not supported for the nmslib search method")
            return self.query_by_sequence_similarity_with_nmslib(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours=max_search_nearest_neighbours, limit_per_sequence=limit_per_sequence)
        elif search_method == 'scann':
            return self.query_by_sequence_similarity_with_scann(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, naive=False, preload_db=preload_db, max_search_nearest_neighbours=max_search_nearest_neighbours, limit_per_sequence=limit_per_sequence)
        elif search_method == 'smafa-naive':
            return self.query_by_sequence_similarity_with_smafa_naive(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, preload_db=preload_db, limit_per_sequence=limit_per_sequence)
        else:
            raise Exception("Unknown search method {}".format(search_method))

    def query_by_sequence_similarity_with_nmslib(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours=None, limit_per_sequence=None):
        logging.info("Searching with nmslib by {} sequence ..".format(sequence_type))

        if max_search_nearest_neighbours is None:
            max_search_nearest_neighbours = max_nearest_neighbours

        last_marker = None
        index = None
        for q in queries:
            if q.marker != last_marker:
                del index
                last_marker = q.marker
                index = sdb.get_sequence_index(last_marker, 'nmslib', sequence_type)
                if index is None:
                    raise Exception("The marker '{}' does not appear to be in the singlem db".format(last_marker))
                logging.info("Querying index for {}".format(last_marker))
                query = select([Marker.id]).where(Marker.marker == last_marker)
                m = sdb.sqlalchemy_connection.execute(query).first()
                if m is None:
                    raise Exception("Marker {} not in the SQL DB".format(last_marker))
                last_marker_id = m['id']

            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                query_protein_sequence = None
                kNN = index.knnQuery(sequence_database.nucleotides_to_binary(q.sequence), max_search_nearest_neighbours)
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                query_protein_sequence = sequence_database.nucleotides_to_protein(q.sequence)
                kNN = index.knnQuery(sequence_database.protein_to_binary(query_protein_sequence), max_search_nearest_neighbours)
            else:
                raise Exception("Unexpected sequence_type")

            num_reported = 0
            for (hit_index, hamming_distance) in zip(kNN[0], kNN[1]):
                div = int(hamming_distance / 2)
                if max_divergence is None or div <= max_divergence:
                    for qres in self.query_result_from_db(sdb, q, sequence_type, hit_index, last_marker, last_marker_id, div, query_protein_sequence=query_protein_sequence, limit_per_sequence=limit_per_sequence):
                        yield qres
                    num_reported += 1
                    if num_reported >= max_nearest_neighbours:
                        break


    def query_by_sequence_similarity_with_scann(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, naive=False, preload_db=False, max_search_nearest_neighbours=None, limit_per_sequence=None):
        if naive:
            logging.info("Searching with SCANN NAIVE by {} sequence ..".format(sequence_type))
        else:
            logging.info("Searching with SCANN by {} sequence ..".format(sequence_type))

        if max_search_nearest_neighbours is None:
            max_search_nearest_neighbours = max_nearest_neighbours

        current_preloaded_db = None

        for marker, marker_queries in itertools.groupby(queries, lambda x: x.marker):
            # Setup index
            if naive:
                index_format = 'scann-naive'
            else:
                index_format = 'scann'
            index = sdb.get_sequence_index(marker, index_format, sequence_type)
            if index is None:
                raise Exception("The marker '{}' does not appear to be '{}' indexed in the singlem db".format(marker, index_format))
            logging.info("Querying index for {}".format(marker))
            query = select([Marker.id]).where(Marker.marker == marker)
            m = sdb.sqlalchemy_connection.execute(query).first()
            if m is None:
                raise Exception("Marker {} not in the SQL DB".format(marker))
            marker_id = m['id']

            # Preload DB if needed
            if preload_db:
                if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                    current_preloaded_db = self.preload_nucleotide_db(sdb, marker_id)
                    current_preloaded_db_indices = pd.Series(
                        current_preloaded_db.groupby('nucleotides_marker_wise_id').indices)
                elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                    current_preloaded_db = self.preload_protein_db(sdb, marker_id)
                    current_preloaded_db_indices = pd.Series(
                        current_preloaded_db.groupby('proteins_marker_wise_id').indices)
                else:
                    raise Exception("Unexpected sequence_type")
                # Convert everything to numpy arrays for speedier indexing than iloc
                current_preloaded_db_sample_name = current_preloaded_db.xs('sample_name',axis=1).to_numpy()
                current_preloaded_db_count = current_preloaded_db.xs('num_hits',axis=1).to_numpy()
                current_preloaded_db_sequence = current_preloaded_db.xs('nucleotide_sequence',axis=1).to_numpy()
                current_preloaded_db_coverage = current_preloaded_db.xs('coverage',axis=1).to_numpy()
                current_preloaded_db_taxonomy = current_preloaded_db.xs('taxonomy',axis=1).to_numpy()
                if sequence_type == SequenceDatabase.PROTEIN_TYPE:
                    current_preloaded_db_protein_sequence = current_preloaded_db.xs('protein_sequence',axis=1).to_numpy()
                del current_preloaded_db
                # current_preloaded_db_count = current_preloaded_db_indices
                if limit_per_sequence:
                    # shuffle and truncate once up front
                    current_preloaded_db_indices.apply(np.random.shuffle)
                    current_preloaded_db_indices = pd.Series([a[:limit_per_sequence] for a in current_preloaded_db_indices])

            # When scann DB is absent due to too few seqs
            if index is None:
                continue

            # Actually do searches, in batches
            for chunked_queries1 in iterable_chunks(marker_queries, 1000):
                chunked_queries = list([a for a in chunked_queries1 if a is not None]) # Remove trailing Nones from the iterable

                if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                    query_array = np.array([
                        np.array(sequence_database.nucleotides_to_binary_array(q.sequence)) for q in chunked_queries])
                elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                    query_protein_sequences = np.array([
                        sequence_database.nucleotides_to_protein(q.sequence) for q in chunked_queries])
                    query_array = np.array([
                        np.array(sequence_database.protein_to_binary_array(seq)) for seq in query_protein_sequences])
                else:
                    raise Exception("Unexpected sequence_type")

                normed = np.array([qa / np.linalg.norm(qa) for qa in query_array])
                kNN_batch = index.search_batched_parallel(normed, max_search_nearest_neighbours)

                for i, q in enumerate(chunked_queries):
                    num_reported = 0
                    for (hit_index, dist) in zip(kNN_batch[0][i], kNN_batch[1][i]):
                        if math.isnan(dist):
                            # Happens when we ask for more hits than there are
                            # in DB in total i.e. for tiny databases.
                            continue

                        if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                            div = round((1.0-float(dist))*len(q.sequence)) # Not sure why this is necessary, why doesn't it return a real distance?
                        else:
                            div = round((1.0-float(dist))*len(query_protein_sequences[0])) # Not sure why this is necessary, why doesn't it return a real distance?

                        if max_divergence is None or div <= max_divergence:
                            if preload_db:
                                hit_index = int(hit_index) # Needed only for tiny databases?
                                if hit_index <= 16 and hit_index >= len(current_preloaded_db_indices):
                                    logging.debug("Skipping hit index {} because it is out of bounds, so a dummy entry".format(hit_index))
                                    continue
                                for entry_i in current_preloaded_db_indices.iat[hit_index]:
                                    otu = OtuTableEntry()
                                    otu.marker = marker
                                    otu.sample_name = current_preloaded_db_sample_name[entry_i]
                                    otu.count = current_preloaded_db_count[entry_i]
                                    otu.sequence = current_preloaded_db_sequence[entry_i]
                                    otu.coverage = current_preloaded_db_coverage[entry_i]
                                    otu.taxonomy = current_preloaded_db_taxonomy[entry_i]
                                    if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                                        yield QueryResult(q, otu, div)
                                    else:
                                        yield QueryResult(
                                            q, otu, div, 
                                            query_protein_sequence=query_protein_sequences[i],
                                            subject_protein_sequence=current_preloaded_db_protein_sequence[entry_i])
                            else:
                                for qres in self.query_result_from_db(sdb, q, sequence_type, hit_index, marker, marker_id, div, 
                                    query_protein_sequence=query_protein_sequences[i] if sequence_type == SequenceDatabase.PROTEIN_TYPE else None,
                                    limit_per_sequence=limit_per_sequence):

                                    yield qres
                            num_reported += 1
                            if num_reported >= max_nearest_neighbours:
                                break

    def query_by_sequence_similarity_with_smafa_naive(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, preload_db=False, max_search_nearest_neighbours=None, limit_per_sequence=None):
        logging.info("Searching with SMAFA NAIVE by {} sequence ..".format(sequence_type))

        if max_nearest_neighbours != 1:
            raise Exception("SMAFA NAIVE only supports max_nearest_neighbours=1 right now")
        if preload_db:
            raise Exception("SMAFA NAIVE does not support preloading the DB right now")

        if max_search_nearest_neighbours is None:
            max_search_nearest_neighbours = max_nearest_neighbours

        current_preloaded_db = None

        for marker, marker_queries in itertools.groupby(queries, lambda x: x.marker):
            index = sdb.get_sequence_index(marker, sequence_database.SMAFA_NAIVE_INDEX_FORMAT, sequence_type)
            if index is None:
                raise Exception("The marker '{}' does not appear to be 'smafa-naive/{}' indexed in the singlem db".format(marker, sequence_type))
            logging.info("Querying index for {}".format(marker))
            query = select([Marker.id]).where(Marker.marker == marker)
            m = sdb.sqlalchemy_connection.execute(query).first()
            if m is None:
                raise Exception("Marker {} not in the SQL DB".format(marker))
            marker_id = m['id']
            
            # Actually do searches, in batches
            for chunked_queries1 in iterable_chunks(marker_queries, 1000):
                chunked_queries = list([a for a in chunked_queries1 if a is not None]) # Remove trailing Nones from the iterable

                if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                    query_array = np.array([
                        np.array(sequence_database.nucleotides_to_binary_array(q.sequence)) for q in chunked_queries])
                elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                    raise NotImplementedError("SMAFA NAIVE does not support protein sequences yet")
                    # query_protein_sequences = np.array([
                    #     sequence_database.nucleotides_to_protein(q.sequence) for q in chunked_queries])
                    # query_array = np.array([
                    #     np.array(sequence_database.protein_to_binary_array(seq)) for seq in query_protein_sequences])
                else:
                    raise Exception("Unexpected sequence_type")

                # Setup a pipe so that the query fasta isn't in a file and can be piped through
                with tempfile.NamedTemporaryFile(prefix='singlem-smafa-naive-stderr-') as f:
                    with subprocess.Popen(
                        ['bash','-c','smafa query --database \'{}\' --query /dev/stdin 2> {}'.format(index, f.name)], stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        bufsize=1,
                        universal_newlines=True) as p:

                        for i, q in enumerate(chunked_queries):
                            num_reported = 0
                            print(">{}\n{}".format(i, q.sequence), file=p.stdin)
                        p.stdin.close()

                        # Read the output
                        previous_query_index = None
                        previous_query_index_num_reported = 0
                        for line in p.stdout:
                            splits = line.strip().split("\t")
                            if len(splits) != 4:
                                raise Exception("Unexpected output from smafa: {}".format(line))
                            query_index_str, hit_index, div_str, _ = splits
                            query_index = int(query_index_str)
                            div = int(div_str)

                            if max_divergence is None or div <= max_divergence:
                                if False: #preload_db:
                                    for entry_i in current_preloaded_db_indices.iat[hit_index]:
                                        otu = OtuTableEntry()
                                        otu.marker = marker
                                        otu.sample_name = current_preloaded_db_sample_name[entry_i]
                                        otu.count = current_preloaded_db_count[entry_i]
                                        otu.sequence = current_preloaded_db_sequence[entry_i]
                                        otu.coverage = current_preloaded_db_coverage[entry_i]
                                        otu.taxonomy = current_preloaded_db_taxonomy[entry_i]
                                        if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                                            yield QueryResult(q, otu, div)
                                        else:
                                            yield QueryResult(
                                                q, otu, div, 
                                                query_protein_sequence=query_protein_sequences[i],
                                                subject_protein_sequence=current_preloaded_db_protein_sequence[entry_i])
                                else:
                                    for qres in self.query_result_from_db(sdb, chunked_queries[query_index], sequence_type, hit_index, marker, marker_id, div, 
                                            query_protein_sequence=query_protein_sequences[i] if sequence_type == SequenceDatabase.PROTEIN_TYPE else None,
                                            limit_per_sequence=limit_per_sequence):

                                            yield qres

                                if query_index == previous_query_index:
                                    previous_query_index_num_reported += 1
                                else:
                                    previous_query_index_num_reported = 0
                                    previous_query_index = query_index
                                if previous_query_index_num_reported > max_nearest_neighbours:
                                    break
                        
                        # Check for errors
                        p.wait()
                        if p.returncode != 0:
                            with open(f.name, 'r') as f:
                                raise Exception("Error running smafa query: {}".format(f.read()))


    def query_by_sequence_similarity_with_annoy(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, max_search_nearest_neighbours=None, limit_per_sequence=None):
        logging.info("Searching with annoy by {} sequence ..".format(sequence_type))

        if max_search_nearest_neighbours is None:
            max_search_nearest_neighbours = max_nearest_neighbours

        last_marker = None
        last_marker_id = None
        index = None
        for q in queries:
            if q.marker != last_marker:
                del index
                last_marker = q.marker
                index = sdb.get_sequence_index(last_marker, 'annoy', sequence_type)
                if index is None:
                    raise Exception("The marker '{}' does not appear to be in the singlem db".format(last_marker))
                logging.info("Querying index for {}".format(last_marker))
                query = select([Marker.id]).where(Marker.marker == last_marker)
                m = sdb.sqlalchemy_connection.execute(query).first()
                if m is None:
                    raise Exception("Marker {} not in the SQL DB".format(last_marker))
                last_marker_id = m['id']

            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                query_protein_sequence = None
                kNN = index.get_nns_by_vector(sequence_database.nucleotides_to_binary_array(q.sequence), max_search_nearest_neighbours, include_distances=True)
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                query_protein_sequence = sequence_database.nucleotides_to_protein(q.sequence)
                kNN = index.get_nns_by_vector(sequence_database.protein_to_binary_array(query_protein_sequence), max_search_nearest_neighbours, include_distances=True)
            else:
                raise Exception("Unexpected sequence_type")

            num_reported = 0
            for (hit_index, hamming_distance) in zip(kNN[0], kNN[1]):
                div = int(hamming_distance / 2)
                if max_divergence is None or div <= max_divergence:
                    for qres in self.query_result_from_db(sdb, q, sequence_type, hit_index, last_marker, last_marker_id, div, query_protein_sequence=query_protein_sequence, limit_per_sequence=limit_per_sequence):
                        yield qres
                    num_reported += 1
                    if num_reported >= max_nearest_neighbours:
                        break

    def query_result_from_db(self, sdb, query, sequence_type, hit_index, marker, marker_id, div, query_protein_sequence=None, limit_per_sequence=None):
        if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
            query2 = select(
                Otu.marker_id, Otu.sample_name, Otu.sequence, Otu.num_hits, Otu.coverage, Otu.taxonomy_id
            ).where(Otu.marker_wise_sequence_id == int(hit_index)).where(Otu.marker_id == int(marker_id)).limit(limit_per_sequence)
            for row in sdb.sqlalchemy_connection.execute(query2):
                otu = OtuTableEntry()
                otu.marker = marker
                otu.sample_name = row.sample_name
                otu.count = row.num_hits
                otu.sequence = row.sequence
                otu.coverage = row.coverage
                otu.taxonomy = sdb.get_taxonomy_via_cache(row.taxonomy_id)
                yield QueryResult(query, otu, div)

        elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
            query2 = select([
                Otu.sequence,
                ProteinSequence.protein_sequence,
                Otu.sample_name,
                Otu.num_hits,
                Otu.coverage,
                Otu.taxonomy_id]) \
                    .where(Otu.taxonomy_id == Taxonomy.id) \
                    .where(Otu.marker_id == marker_id) \
                    .where(NucleotidesProteins.nucleotide_id == Otu.sequence_id) \
                    .where(NucleotidesProteins.protein_id == ProteinSequence.id) \
                    .where(ProteinSequence.marker_wise_id == int(hit_index)) \
                    .where(Otu.marker_id == marker_id)
            if limit_per_sequence is not None:
                query2 = query2.limit(limit_per_sequence)

            results = sdb.sqlalchemy_connection.execute(query2)

            if results is None and hit_index <= 16:
                # For very small indexes, SCANN can have dummy sequences that
                # are not in the SQL DB. Ignore these.
                pass
            else:
                for entry in results:
                    otu = OtuTableEntry()
                    otu.marker = marker
                    otu.sample_name = entry['sample_name']
                    otu.sequence = entry['sequence']
                    otu.count = entry['num_hits']
                    otu.coverage = entry['coverage']
                    otu.taxonomy = sdb.get_taxonomy_via_cache(entry['taxonomy_id'])
                    yield QueryResult(query, otu, div, query_protein_sequence=query_protein_sequence, subject_protein_sequence=entry['protein_sequence'])
        else:
            raise Exception("unknown sequence_type")

    def divergence(self, seq1, seq2):
        """Return the number of bases two sequences differ by"""
        if len(seq1) != len(seq2):
            raise Exception(
                "Attempted comparison of two OTU sequences with differing lengths: {} and {}" \
                .format(seq1, seq2))
        d = 0
        for (a, b) in zip(seq1, seq2):
            if a != b: d += 1
        return d


    def query_by_sqlite(self, queries, db):
        max_set_size = 999 # Cannot query sqlite with > 999 '?' entries, so
                           # query in batches.
        sequence_to_query_id = {}
        queries_list = list(queries)
        seqs = set()
        for i, query in enumerate(queries_list):
            seqs.add(query.sequence)
            try:
                sequence_to_query_id[query.sequence].append(i)
            except KeyError:
                sequence_to_query_id[query.sequence] = [i]

        for chunk in iterable_chunks(seqs, max_set_size):
            with db.engine.connect() as connection:
                stmt = select([
                    Otu.sample_name,
                    Otu.num_hits,
                    Otu.coverage,
                    Otu.taxonomy_id,
                    NucleotideSequence.sequence,
                    Marker.marker,
                    ]).select_from(Otu).join(NucleotideSequence).join(Marker).filter(NucleotideSequence.sequence.in_([seq for seq in chunk if seq is not None]))
                for entry in connection.execute(stmt):
                    for qid in sequence_to_query_id[entry.sequence]:
                        otu = OtuTableEntry()
                        otu.marker = entry.marker
                        otu.sample_name = entry.sample_name
                        otu.sequence = entry.sequence
                        otu.count = entry.num_hits
                        otu.coverage = entry.coverage
                        otu.taxonomy = db.get_taxonomy_via_cache(entry.taxonomy_id)
                        yield QueryResult(queries_list[qid], otu, 0)

    def print_samples(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'), min_version=5)
        sample_names = kwargs.pop('sample_names')
        taxonomy = kwargs.pop('taxonomy')
        output_io = kwargs.pop('output_io')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        max_set_size = 999 # Cannot query sqlite with > 999 '?' entries, so
                           # query in batches.
        if sample_names:
            query_chunks = set(sample_names)
        else:
            query_chunks = [taxonomy]
        otus = OtuTable()
        total_printed = 0
        first_chunk = True
        db.engine.execution_options(stream_results=True) # So RAM usage is not crazy for large queries
        for chunk in SequenceDatabase._grouper(query_chunks, max_set_size):
            if sample_names:
                query = select(Otu).filter(Otu.sample_name.in_([sample for sample in chunk if sample is not None]))
            elif taxonomy:
                query = select(Otu) \
                    .where(Taxonomy.taxonomy.like('%{}%'.format(taxonomy))) \
                    .where(Otu.taxonomy_id == Taxonomy.id)
            else:
                raise Exception("Programming error")

            for entry in db.sqlalchemy_connection.execute(query):
                otu = OtuTableEntry()
                otu.marker = db.get_marker_via_cache(entry['marker_id'])
                otu.sample_name = entry['sample_name']
                otu.sequence = entry['sequence']
                otu.count = entry['num_hits']
                otu.coverage = entry['coverage']
                otu.taxonomy = db.get_taxonomy_via_cache(entry['taxonomy_id'])
                otus.add([otu])
                total_printed += 1
            otus.write_to(output_io, print_header=first_chunk)
            otus = OtuTable()
            first_chunk = False
        logging.info("Printed %i OTU table entries" % total_printed)


class QueryInputSequence:
    def __init__(self, name, sequence, marker):
        self.name = name
        self.sequence = sequence
        self.marker = marker

class QueryResult:
    def __init__(self, query, subject, divergence, query_protein_sequence=None, subject_protein_sequence=None):
        self.query = query
        self.subject = subject
        self.divergence = divergence
        self.query_protein_sequence = query_protein_sequence
        self.subject_protein_sequence = subject_protein_sequence

class MarkerSortedQueryInput:
    def __init__(self, sorted_io):
        self._original_io = sorted_io
        self._reopened_io = open(sorted_io.name)
        self._reopened_reader = iter(csv.reader(self._reopened_io, delimiter="\t"))

    def __iter__(self):
        return self

    def __next__(self):
        '''implement iterable'''
        row = next(self._reopened_reader) #propogate the StopIteration
        if len(row) != 3:
            raise Exception("Strange query formating issue")
        return QueryInputSequence(row[2], row[1], row[0])

    def __del__(self):
        self._original_io.close()
        self._reopened_io.close()
