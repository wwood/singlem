import tempfile
import logging
import subprocess
import sys
import os
from orator import DatabaseManager, Model
from orator.exceptions.query import QueryException
from Bio import pairwise2
import numpy as np
import csv

from .sequence_database import SequenceDatabase
from . import sequence_database
from .sequence_classes import SeqReader
from .query_formatters import SparseResultFormatter
from .otu_table_collection import OtuTableCollection
from .otu_table_entry import OtuTableEntry
from .otu_table import OtuTable
from .otu_table_collection import StreamingOtuTableCollection

class Querier:
    def query(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'))
        max_divergence = kwargs.pop('max_divergence')
        output_style = kwargs.pop('output_style')
        query_otu_table = kwargs.pop('query_otu_table')
        num_threads = kwargs.pop('num_threads')
        search_method = kwargs.pop('search_method')
        sequence_type = kwargs.pop('sequence_type')
        # stream_output = kwargs.pop('stream_output')
        max_nearest_neighbours = kwargs.pop('max_nearest_neighbours')
        
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        queries = self.prepare_query_sequences(query_otu_table, num_threads)

        query_results = self.query_with_queries(
            queries, db, max_divergence, search_method, sequence_type, max_nearest_neighbours=max_nearest_neighbours)
        # Only reason not to stream would be so that the queries are returned in the same order as passed in. eh for now.
        do_stream = True
        if not do_stream:
            query_results = list(query_results)
            logging.info("Printing %i hits" % len(query_results))

        if output_style == 'sparse':
            SparseResultFormatter().write(query_results, sys.stdout, streaming=do_stream)
        elif output_style == None:
            return query_results
        else:
            raise Exception("Programming error")

    def query_subject_otu_table(self, **kwargs):
        subject_otus = kwargs.pop('subject_otu_collection')
        query_sequence = kwargs.pop('query_sequence')
        max_divergence = kwargs.pop('max_divergence')
        query_otu_table = kwargs.pop('query_otu_table')
        query_fasta = kwargs.pop('query_fasta')
        output_style = kwargs.pop('output_style')

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if (query_otu_table and query_sequence) or \
            (query_otu_table and query_fasta) or \
            (query_sequence and query_fasta):
            raise Exception("Only one of --query_fasta, --query_otu_table and --query_sequence is allowable")

        queries = self.prepare_query_sequences(
            query_sequence, query_otu_table, query_fasta)
        logging.info("Read in %i queries" % len(queries))

        query_results = []
        for query in queries:
            for otu in subject_otus:
                logging.debug("Comparing query {}/{} with subject {}/{}".format(
                    query.name,
                    query.sequence,
                    otu.sample_name,
                    otu.sequence))
                # TODO: Optimise these alignment scoring parameters for a codon-wise alignment?
                alignments = pairwise2.align.globalxs(
                    query.sequence.upper(),
                    otu.sequence.upper(),
                    0, 0) # Use 0 not -1 because the gap is counted on top of the mismatch
                best_alignment = max(alignments, key=lambda p: p[2])
                divergence = max(len(query.sequence), len(otu.sequence)) - int(best_alignment[2])
                logging.debug("Found divergence of best alignment {}".format(divergence))
                if divergence <= max_divergence:
                    query_results.append(QueryResult(query, otu, divergence))

        if output_style == 'sparse':
            SparseResultFormatter().write(query_results, sys.stdout)
        elif output_style == None:
            return query_results
        else:
            raise Exception("Programming error")



    def prepare_query_sequences(self, query_otu_table, num_threads):
        '''return an iterable of QueryInputSequence objects sorted by marker.'''
        otus = StreamingOtuTableCollection()
        otus.add_otu_table_file(query_otu_table)

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

    def query_with_queries(self, queries, sdb, max_divergence, search_method, sequence_type, max_nearest_neighbours):
        sdb.query_builder(check=True)
        if max_divergence == 0 and sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
            return self.query_by_sqlite(queries, sdb)
        elif search_method == 'naive':
            return self.query_by_sequence_similarity_with_scann(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, naive=True)
        elif search_method == 'annoy':
            return self.query_by_sequence_similarity_with_annoy(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours)
        elif search_method == 'nmslib':
            return self.query_by_sequence_similarity_with_nmslib(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours)
        elif search_method == 'scann':
            return self.query_by_sequence_similarity_with_scann(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, naive=False)
        else:
            raise Exception("Unknown search method {}".format(search_method))

    ### Method no longer used, scann is used instead
    # def query_by_sequence_similarity_naively(self, queries, sdb, max_divergence, sequence_type):
    #     if sequence_type == SequenceDatabase.PROTEIN_TYPE:
    #         raise Exception("Naive search by amino acid not implemented for now")
    #     marker_to_queries = {}
    #     results = []
    #     logging.info("Searching by nucleotide sequence with (unoptimised) naive search ..")

    #     for query in queries:
    #         if query.marker not in marker_to_queries:
    #             marker_to_queries[query.marker] = []
    #         marker_to_queries[query.marker].append(query)

    #     for marker, subqueries in marker_to_queries.items():
    #         for (i, entry) in enumerate(sdb.query_builder().table('otus'). \
    #             join('markers','marker_id','=','markers.id').where('markers.marker',marker). \
    #             join('nucleotides','sequence_id','=','nucleotides.id'). \
    #             get()):
    #             if len(subqueries) >= 5 and i > 0 and i % 100000 == 0:
    #                 logging.info("Processing query {} for marker {}".format(i, marker))

    #             for q in subqueries:
    #                 div = self.divergence(q.sequence, entry['sequence'])
    #                 if div <= max_divergence:
    #                     otu = OtuTableEntry()
    #                     otu.marker = entry['marker']
    #                     otu.sample_name = entry['sample_name']
    #                     otu.sequence = entry['sequence']
    #                     otu.count = entry['num_hits']
    #                     otu.coverage = entry['coverage']
    #                     otu.taxonomy = entry['taxonomy']
    #                     yield QueryResult(q, otu, div)


    def query_by_sequence_similarity_with_nmslib(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours):
        marker_to_queries = {}
        results = []
        logging.info("Searching with nmslib by {} sequence ..".format(sequence_type))

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
                
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                kNN = index.knnQuery(sequence_database.nucleotides_to_binary(q.sequence), max_nearest_neighbours)
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                kNN = index.knnQuery(sequence_database.protein_to_binary(sequence_database.nucleotides_to_protein(q.sequence)), max_nearest_neighbours)
            else:
                raise Exception("Unexpected sequence_type")

            for (hit_index, hamming_distance) in zip(kNN[0], kNN[1]):
                div = int(hamming_distance / 2)
                if div <= max_divergence:
                    for qres in self.query_result_from_db(sdb, q, sequence_type, hit_index, last_marker, div):
                        yield qres


    def query_by_sequence_similarity_with_scann(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours, naive=False):
        logging.info("Searching with SCANN by {} sequence ..".format(sequence_type))

        last_marker = None
        last_marker_id = None
        index = None
        for q in queries:
            if q.marker != last_marker:
                del index
                last_marker = q.marker
                if naive:
                    index_format = 'scann-naive'
                else:
                    index_format = 'scann'
                index = sdb.get_sequence_index(last_marker, index_format, sequence_type)
                logging.info("Querying index for {}".format(last_marker))
                m = sdb.query_builder().table('markers').where('marker',last_marker).first()
                if m is None:
                    raise Exception("Marker {} not in the SQL DB".format(last_marker))
                last_marker_id = m['id']
            # When scann DB is absent due to too few seqs
            if index is None:
                continue

            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                query_array = np.array(sequence_database.nucleotides_to_binary_array(q.sequence))
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                query_array = np.array(
                    sequence_database.protein_to_binary_array(
                        sequence_database.nucleotides_to_protein(q.sequence)))
            else:
                raise Exception("Unexpected sequence_type")
            normed = query_array / np.linalg.norm(query_array)
            kNN = index.search(normed, max_nearest_neighbours)

            for (hit_index, dist) in zip(kNN[0], kNN[1]): # Possibly can know div distance from scann distance so less SQL?
                div = int((1.0-dist)*len(q.sequence)) # Not sure why this is necessary, why doesn't it return a real distance?
                if div <= max_divergence:
                    for qres in self.query_result_from_db(sdb, q, sequence_type, hit_index, last_marker, div):
                        yield qres

    def query_result_from_db(self, sdb, query, sequence_type, hit_index, marker, div):
        # TODO: Could be faster by caching the marker_id?
        if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
            for entry in sdb.query_builder().table('otus'). \
                join('markers','marker_id','=','markers.id').where('markers.marker',marker). \
                join('nucleotides','sequence_id','=','nucleotides.id'). \
                where('nucleotides.marker_wise_id', int(hit_index)). \
                get():

                otu = OtuTableEntry()
                otu.marker = entry['marker']
                otu.sample_name = entry['sample_name']
                otu.sequence = entry['sequence']
                otu.count = entry['num_hits']
                otu.coverage = entry['coverage']
                otu.taxonomy = entry['taxonomy']
                yield QueryResult(query, otu, div)
        elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
            for entry in sdb.query_builder().table('otus'). \
                join('markers','marker_id','=','markers.id').where('markers.marker',marker). \
                join('nucleotides','sequence_id','=','nucleotides.id'). \
                join('nucleotides_proteins','nucleotides_proteins.nucleotide_id','=','nucleotides.id'). \
                join('proteins','nucleotides_proteins.protein_id','=','proteins.id'). \
                where('proteins.marker_wise_id', int(hit_index)). \
                get():

                otu = OtuTableEntry()
                otu.marker = entry['marker']
                otu.sample_name = entry['sample_name']
                otu.sequence = entry['sequence']
                otu.count = entry['num_hits']
                otu.coverage = entry['coverage']
                otu.taxonomy = entry['taxonomy']
                yield QueryResult(query, otu, div)
        else:
            raise Exception("unknown sequence_type")

        

    def query_by_sequence_similarity_with_annoy(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours):
        logging.info("Searching with annoy by {} sequence ..".format(sequence_type))

        last_marker = None
        index = None
        for q in queries:
            if q.marker != last_marker:
                del index
                last_marker = q.marker
                index = sdb.get_sequence_index(last_marker, 'annoy', sequence_type)
                if index is None:
                    raise Exception("The marker '{}' does not appear to be in the singlem db".format(last_marker))
                logging.info("Querying index for {}".format(last_marker))

            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                kNN = index.get_nns_by_vector(sequence_database.nucleotides_to_binary_array(q.sequence), max_nearest_neighbours, include_distances=True)
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                kNN = index.get_nns_by_vector(sequence_database.protein_to_binary_array(sequence_database.nucleotides_to_protein(q.sequence)), max_nearest_neighbours, include_distances=True)
            else:
                raise Exception("Unexpected sequence_type")

            for (hit_index, hamming_distance) in zip(kNN[0], kNN[1]):
                div = int(hamming_distance / 2)
                if div <= max_divergence:
                    for qres in self.query_result_from_db(sdb, q, sequence_type, hit_index, last_marker, div):
                        yield qres

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

        results = []
        for chunk in SequenceDatabase.grouper(seqs, max_set_size):
            for entry in db.table('otus').where_in(
                    'sequence', [seq for seq in chunk if seq is not None]).get():
                for qid in sequence_to_query_id[entry.sequence]:
                    otu = OtuTableEntry()
                    otu.marker = entry.marker
                    otu.sample_name = entry.sample_name
                    otu.sequence = entry.sequence
                    otu.count = entry.num_hits
                    otu.coverage = entry.coverage
                    otu.taxonomy = entry.taxonomy
                    yield QueryResult(queries_list[qid], otu, 0)

    def print_samples(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'))
        sample_names = kwargs.pop('sample_names')
        taxonomy = kwargs.pop('taxonomy')
        output_io = kwargs.pop('output_io')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        dbm = db.query_builder()

        max_set_size = 999 # Cannot query sqlite with > 999 '?' entries, so
                           # query in batches.
        if sample_names:
            query_chunks = set(sample_names)
        else:
            query_chunks = [taxonomy]
        otus = OtuTable()
        total_printed = 0
        for chunk in SequenceDatabase._grouper(query_chunks, max_set_size):
            if sample_names:
                it = dbm.table('otus').join('markers','marker_id','=','markers.id').join('nucleotides','sequence_id','=','nucleotides.id').where_in(
                    'sample_name', [sample for sample in chunk if sample is not None]).get()
            elif taxonomy:
                it = dbm.table('otus').join('markers','marker_id','=','markers.id').join('nucleotides','sequence_id','=','nucleotides.id').where(
                    'taxonomy', 'like', "%%%s%%" % taxonomy).get()
            else:
                raise Exception("Programming error")

            for entry in it:
                otu = OtuTableEntry()
                otu.marker = entry['marker']
                otu.sample_name = entry['sample_name']
                otu.sequence = entry['sequence']
                otu.count = entry['num_hits']
                otu.coverage = entry['coverage']
                otu.taxonomy = entry['taxonomy']
                otus.add([otu])
                total_printed += 1
        otus.write_to(output_io)
        logging.info("Printed %i OTU table entries" % total_printed)


class QueryInputSequence:
    def __init__(self, name, sequence, marker):
        self.name = name
        self.sequence = sequence
        self.marker = marker

class QueryResult:
    def __init__(self, query, subject, divergence):
        self.query = query
        self.subject = subject
        self.divergence = divergence

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