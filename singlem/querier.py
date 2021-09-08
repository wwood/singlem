import tempfile
import logging
import subprocess
import sys
import os
from orator import DatabaseManager, Model
from orator.exceptions.query import QueryException
from Bio import pairwise2

from .sequence_database import SequenceDatabase
from . import sequence_database
from .sequence_classes import SeqReader
from .query_formatters import SparseResultFormatter
from .otu_table_collection import OtuTableCollection
from .otu_table_entry import OtuTableEntry
from .otu_table import OtuTable

class Querier:
    def query(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'))
        query_sequence = kwargs.pop('query_sequence')
        max_divergence = kwargs.pop('max_divergence')
        output_style = kwargs.pop('output_style')
        query_otu_table = kwargs.pop('query_otu_table')
        query_fasta = kwargs.pop('query_fasta')
        num_threads = kwargs.pop('num_threads')
        search_method = kwargs.pop('search_method')
        sequence_type = kwargs.pop('sequence_type')
        stream_output = kwargs.pop('stream_output')
        max_nearest_neighbours = kwargs.pop('max_nearest_neighbours')
        
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if (query_otu_table and query_sequence) or \
            (query_otu_table and query_fasta) or \
            (query_sequence and query_fasta):
            raise Exception("Only one of --query_fasta, --query_otu_table and --query_sequence is allowable")

        queries = self.prepare_query_sequences(
            query_sequence, query_otu_table, query_fasta)
        logging.info("Read in %i queries" % len(queries))

        query_results = self.query_with_queries(
            queries, db, max_divergence, search_method, sequence_type, max_nearest_neighbours=max_nearest_neighbours)
        # nmslib and annoy are in order per-query already, so we may as well stream
        do_stream = stream_output or search_method in ['nmslib','annoy']
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



    def prepare_query_sequences(self, query_sequence, query_otu_table, query_fasta):
        '''Given potential ways to define query sequences (as file path strings),
        return a list of QueryInputSequence objects.

        '''
        if query_sequence:
            queries = [QueryInputSequence('unnamed_sequence',query_sequence)]
        elif query_otu_table:
            queries = []
            otus = OtuTableCollection()
            with open(query_otu_table) as f:
                otus.add_otu_table(f)
                for e in otus:
                    queries.append(QueryInputSequence(
                        ';'.join([e.sample_name]),
                        e.sequence,
                        e.marker))
        elif query_fasta:
            queries = []
            with open(query_fasta) as f:
                for name, seq, _ in SeqReader().readfq(f):
                    queries.append(QueryInputSequence(
                        name, seq))
        else:
            raise Exception("No query option specified, cannot continue")
        return queries

    def query_with_queries(self, queries, sdb, max_divergence, search_method, sequence_type, max_nearest_neighbours):
        sdb.query_builder(check=True)
        if max_divergence == 0 and sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
            return self.query_by_sqlite(queries, sdb)
        elif search_method == 'naive':
            return self.query_by_sequence_similarity_naively(
                queries, sdb, max_divergence, sequence_type)
        elif search_method == 'annoy':
            return self.query_by_sequence_similarity_with_annoy(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours)
        elif search_method == 'nmslib':
            return self.query_by_sequence_similarity_with_nmslib(
                queries, sdb, max_divergence, sequence_type, max_nearest_neighbours)
        else:
            raise Exception("Unknown search method {}".format(search_method))

    def query_by_sequence_similarity_naively(self, queries, sdb, max_divergence, sequence_type):
        if sequence_type == SequenceDatabase.PROTEIN_TYPE:
            raise Exception("Naive search by amino acid not implemented for now")
        marker_to_queries = {}
        results = []
        logging.info("Searching by nucleotide sequence with (unoptimised) naive search ..")

        for query in queries:
            if query.marker not in marker_to_queries:
                marker_to_queries[query.marker] = []
            marker_to_queries[query.marker].append(query)

        for marker, subqueries in marker_to_queries.items():
            for (i, entry) in enumerate(sdb.query_builder().table('otus'). \
                join('markers','marker_id','=','markers.id').where('markers.marker',marker). \
                join('nucleotides','sequence_id','=','nucleotides.id'). \
                get()):
                if len(subqueries) >= 5 and i > 0 and i % 100000 == 0:
                    logging.info("Processing query {} for marker {}".format(i, marker))

                for q in subqueries:
                    div = self.divergence(q.sequence, entry['sequence'])
                    if div <= max_divergence:
                        otu = OtuTableEntry()
                        otu.marker = entry['marker']
                        otu.sample_name = entry['sample_name']
                        otu.sequence = entry['sequence']
                        otu.count = entry['num_hits']
                        otu.coverage = entry['coverage']
                        otu.taxonomy = entry['taxonomy']
                        yield QueryResult(query, otu, div)


    def query_by_sequence_similarity_with_nmslib(self, queries, sdb, max_divergence, sequence_type, max_nearest_neighbours):
        marker_to_queries = {}
        results = []
        logging.info("Searching with nmslib by nucleotide sequence ..")

        for query in queries:
            if query.marker not in marker_to_queries:
                marker_to_queries[query.marker] = []
            marker_to_queries[query.marker].append(query)

        for marker, subqueries in marker_to_queries.items():
            index = sdb.get_sequence_index(marker, 'nmslib', sequence_type)
            if index is None:
                raise Exception("The marker '{}' does not appear to be in the singlem db".format(marker))
            logging.info("Querying index for {}".format(marker))

            for q in subqueries:
                if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                    kNN = index.knnQuery(sequence_database.nucleotides_to_binary(q.sequence), max_nearest_neighbours)
                elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                    kNN = index.knnQuery(sequence_database.protein_to_binary(sequence_database.nucleotides_to_protein(q.sequence)), max_nearest_neighbours)
                else:
                    raise Exception("Unexpected sequence_type")

                for (hit_index, hamming_distance) in zip(kNN[0], kNN[1]):
                    div = int(hamming_distance / 2)
                    if div <= max_divergence:
                        for qres in self.query_result_from_db(sdb, query, sequence_type, hit_index, marker, div):
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
        marker_to_queries = {}
        results = []
        logging.info("Searching with annoy by {} sequence ..".format(sequence_type))

        for query in queries:
            if query.marker not in marker_to_queries:
                marker_to_queries[query.marker] = []
            marker_to_queries[query.marker].append(query)

        for marker, subqueries in marker_to_queries.items():
            index = sdb.get_sequence_index(marker, 'annoy', sequence_type)
            if index is None:
                raise Exception("The marker '{}' does not appear to be in the singlem db".format(marker))
            logging.info("Querying index for {}".format(marker))

            for q in subqueries:
                if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                    kNN = index.get_nns_by_vector(sequence_database.nucleotides_to_binary_array(q.sequence), max_nearest_neighbours, include_distances=True)
                elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                    kNN = index.get_nns_by_vector(sequence_database.protein_to_binary_array(sequence_database.nucleotides_to_protein(q.sequence)), max_nearest_neighbours, include_distances=True)
                else:
                    raise Exception("Unexpected sequence_type")

                for (hit_index, hamming_distance) in zip(kNN[0], kNN[1]):
                    div = int(hamming_distance / 2)
                    if div <= max_divergence:
                        for qres in self.query_result_from_db(sdb, query, sequence_type, hit_index, marker, div):
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
