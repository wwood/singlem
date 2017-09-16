import tempfile
import extern
import logging
import subprocess
import sys
import os
from orator import DatabaseManager, Model

from sequence_database import SequenceDatabase
from sequence_classes import SeqReader
from query_formatters import SparseResultFormatter
from otu_table_collection import OtuTableCollection
from otu_table_entry import OtuTableEntry

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

        query_results = self.query_with_queries(queries, db, max_divergence)

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


    def query_with_queries(self, queries, db, max_divergence):
        sqlite_db_path = db.sqlite_file
        if not os.path.exists(sqlite_db_path):
            raise Exception("Sqlite database not found at '%s', indicating that either the SingleM database was built with an out-dated SingleM version, or that the database is corrupt. Please generate a new database with the current version of SingleM.")
        logging.info("Connecting to %s" % sqlite_db_path)
        dbm = DatabaseManager({
        'sqlite3': {
            'driver': 'sqlite',
            'database': sqlite_db_path
        }})
        Model.set_connection_resolver(dbm)

        if max_divergence == 0:
            return self.query_by_sqlite(queries, dbm)
        else:
            return self.query_by_smafa(
                queries, db.smafa_dbs(), dbm, max_divergence)


    def query_by_smafa(self, queries, smafa_dbs, sqlite_db,  max_divergence):
        # Generate a tempfile of all the queries
        with tempfile.NamedTemporaryFile(prefix='singlem_query_smafa') as infile:
            sequence_to_queries = {}
            for query in queries:
                infile.write(">%s\n" % query.name)
                infile.write(query.sequence+"\n")
                infile.flush()
                if query.sequence in sequence_to_queries:
                    sequence_to_queries[query.sequence].append(query)
                else:
                    sequence_to_queries[query.sequence] = [query]

            results = []
            for smafa_db in smafa_dbs:
                cmd = "smafa query -q -d %i '%s' '%s'" % (max_divergence, smafa_db, infile.name)
                logging.debug("Running cmd with popen: %s" % cmd)
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                for line in iter(proc.stdout.readline,''):
                    query_name, query_sequence, subject_sequence, divergence = line.split("\t")[:4]
                    queries = sequence_to_queries[query_sequence]
                    for entry in sqlite_db.table('otus').where('sequence',subject_sequence).get():
                        otu = OtuTableEntry()
                        otu.marker = entry.marker
                        otu.sample_name = entry.sample_name
                        otu.sequence = entry.sequence
                        otu.count = entry.num_hits
                        otu.coverage = entry.coverage
                        otu.taxonomy = entry.taxonomy
                        for query in queries:
                            results.append(QueryResult(query, otu, divergence))
        return results


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


class QueryInputSequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

class BlastQueryResultLine:
    def __init__(self, blast_output_line):
        self.qseqid, self.sseqid, _, _, self.mismatch, self.gaps, self.qstart,\
            self.qend, self.sstart, \
            self.send = blast_output_line.strip().split("\t")[:10]

class QueryResult:
    def __init__(self, query, subject, divergence):
        self.query = query
        self.subject = subject
        self.divergence = divergence

class BlastQueryResults:
    def __init__(self, queries_subjects_divergences, sseqid_to_hits):
        self._queries_subjects_divergences = queries_subjects_divergences
        self._sseqid_to_hits = sseqid_to_hits

    def __iter__(self):
        '''Iterate over results, yielding QueryResult objects'''
        for qsd in self._queries_subjects_divergences:
            for subject in self._sseqid_to_hits[qsd[1]]:
                yield QueryResult(qsd[0], subject, qsd[2])
