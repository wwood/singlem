import tempfile
import logging
import subprocess
import sys
import os
from orator import DatabaseManager, Model
from orator.exceptions.query import QueryException
from Bio import pairwise2

from .sequence_database import SequenceDatabase
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
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if (query_otu_table and query_sequence) or \
            (query_otu_table and query_fasta) or \
            (query_sequence and query_fasta):
            raise Exception("Only one of --query_fasta, --query_otu_table and --query_sequence is allowable")

        queries = self.prepare_query_sequences(
            query_sequence, query_otu_table, query_fasta)
        logging.info("Read in %i queries" % len(queries))

        query_results = self.query_with_queries(queries, db, max_divergence)
        logging.info("Printing %i hits" % len(query_results))

        if output_style == 'sparse':
            SparseResultFormatter().write(query_results, sys.stdout)
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
                        ';'.join([e.sample_name,e.marker]),
                        e.sequence))
        elif query_fasta:
            queries = []
            with open(query_fasta) as f:
                for name, seq, _ in SeqReader().readfq(f):
                    queries.append(QueryInputSequence(
                        name, seq))
        else:
            raise Exception("No query option specified, cannot continue")
        return queries

    def _connect_to_sqlite(self, db):
        sqlite_db_path = db.sqlite_file
        if not os.path.exists(sqlite_db_path):
            raise Exception("Sqlite database not found at '%s', indicating that either the SingleM database was built with an out-dated SingleM version, or that the database is corrupt. Please generate a new database with the current version of SingleM.")
        logging.debug("Connecting to %s" % sqlite_db_path)
        dbm = DatabaseManager({
        'sqlite3': {
            'driver': 'sqlite',
            'database': sqlite_db_path
        }})
        Model.set_connection_resolver(dbm)
        try:
            len(dbm.table('otus').limit(1).get())
        except Exception as e:
            logging.error("Failure to extract any data from the otus table of the SQ Lite DB indicates this SingleM DB is either too old or is corrupt.")
            raise(e)
        try:
            len(dbm.table('clusters').limit(1).get())
        except QueryException:
            logging.error("Failure to extract any data from the 'clusters' table indicates this SingleM DB is out-dated, and cannot be used with query implemented in this version of SingleM")
            sys.exit(1)
        return dbm

    def query_with_queries(self, queries, db, max_divergence):
        dbm = self._connect_to_sqlite(db)
        if max_divergence == 0:
            return self.query_by_sqlite(queries, dbm)
        else:
            return self.query_by_smafa(
                queries, db, dbm, max_divergence)


    def query_by_smafa(self, queries, sdb, sqlite_db,  max_divergence):
        # Generate a tempfile of all the queries
        with tempfile.NamedTemporaryFile(prefix='singlem_query_smafa',mode='w') as infile:
            for query in queries:
                infile.write(">%s\n" % query.name)
                infile.write(query.sequence+"\n")
                infile.flush()
            infile.flush()

            results = []
            for smafa_db in sdb.smafa_dbs():
                logging.info("Querying smafadb {}".format(smafa_db))
                # Query with at least the divergence of the clustering + max_divergence,
                # otherwise there may be false negatives.
                smafa_divergence = max_divergence
                min_divergence = sdb.smafa_clustering_divergence() + max_divergence
                if smafa_divergence < min_divergence:
                    smafa_divergence = min_divergence
                cmd = "smafa query -q -d %i '%s' '%s'" % (smafa_divergence, smafa_db, infile.name)
                logging.debug("Running cmd with popen: %s" % cmd)
                proc = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    universal_newlines=True)
                for line in proc.stdout:
                    query_name, query_sequence, subject_sequence, divergence = line.split("\t")[:4]
                    logging.debug("smafa query with {} returned hit sequence {}".format(
                        query_name, subject_sequence))
                    # Query for OTUs that have a cluster representative that is the hit sequence
                    for entry in sqlite_db \
                        .table('otus') \
                        .join('clusters','clusters.member','=','otus.sequence') \
                        .where('clusters.representative',subject_sequence).get():
                        logging.debug(
                            "Returned a potential hit through through clustering table: {} (divergence {})".format(
                                entry.sequence, divergence))

                        if sdb.smafa_clustering_divergence() == 0:
                            div = int(divergence)
                        else:
                            div = self.divergence(query_sequence, entry.sequence)
                            logging.debug("Calculated divergence to be {}".format(div))
                        if div <= max_divergence:
                            #NOTE: Not the same objects as the original query objects.
                            query = QueryInputSequence(query_name, query_sequence)
                            otu = OtuTableEntry()
                            otu.marker = entry.marker
                            otu.sample_name = entry.sample_name
                            otu.sequence = entry.sequence
                            otu.count = entry.num_hits
                            otu.coverage = entry.coverage
                            otu.taxonomy = entry.taxonomy
                            results.append(QueryResult(query, otu, div))
        return results

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
                    results.append(QueryResult(queries_list[qid], otu, 0))
        return results

    def print_samples(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'))
        sample_names = kwargs.pop('sample_names')
        taxonomy = kwargs.pop('taxonomy')
        output_io = kwargs.pop('output_io')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        dbm = self._connect_to_sqlite(db)

        max_set_size = 999 # Cannot query sqlite with > 999 '?' entries, so
                           # query in batches.
        if sample_names:
            query_chunks = set(sample_names)
        else:
            query_chunks = [taxonomy]
        otus = OtuTable()
        total_printed = 0
        for chunk in SequenceDatabase.grouper(query_chunks, max_set_size):
            if sample_names:
                it = dbm.table('otus').where_in(
                    'sample_name', [sample for sample in chunk if sample is not None]).get()
            elif taxonomy:
                it = dbm.table('otus').where(
                    'taxonomy', 'like', "%%%s%%" % taxonomy).get()
            else:
                raise Exception("Programming error")

            for entry in it:
                otu = OtuTableEntry()
                otu.marker = entry.marker
                otu.sample_name = entry.sample_name
                otu.sequence = entry.sequence
                otu.count = entry.num_hits
                otu.coverage = entry.coverage
                otu.taxonomy = entry.taxonomy
                otus.add([otu])
                total_printed += 1
        otus.write_to(output_io)
        logging.info("Printed %i OTU table entries" % total_printed)


class QueryInputSequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

class QueryResult:
    def __init__(self, query, subject, divergence):
        self.query = query
        self.subject = subject
        self.divergence = divergence
