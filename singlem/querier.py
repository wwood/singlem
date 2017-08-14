import tempfile
import extern
import logging
import subprocess
import sys
import os
from orator import DatabaseManager, Model

from sequence_database import SequenceDatabase
from sequence_classes import SeqReader
from query_formatters import SparseResultFormatter#, DenseResultFormatter
from otu_table_collection import OtuTableCollection

class Querier:
    def query(self, **kwargs):
        db = SequenceDatabase.acquire(kwargs.pop('db'))
        query_sequence = kwargs.pop('query_sequence')
        max_target_seqs = kwargs.pop('max_target_seqs')
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

        if max_divergence == 0:
            sqlite_db_path = db.sqlite_file
            if not os.path.exists(sqlite_db_path):
                raise Exception("Sqlite database not found at '%s', indicating that either the SingleM database was built with an out-dated SingleM version, or that the database is corrupt. Please generate a new database with the current version of SingleM.")
            query_results = self.query_by_sqlite(queries, sqlite_db_path)
        else:
            query_results = self.query_by_blast(
                queries, db, max_target_seqs, max_divergence, num_threads)

        if output_style == 'sparse':
            SparseResultFormatter().write(query_results, sys.stdout)
        elif output_style == None:
            return query_results
        else:
            raise Exception()


    def query_by_blast(self, queries, db, max_target_seqs, max_divergence, num_threads):
        # blast the query against the database, output as csv
        found_distances_and_names = []
        with tempfile.NamedTemporaryFile(prefix='singlem_query') as infile:
            for i, query in enumerate(queries):
                infile.write(">%i\n" % i)
                infile.write(query.sequence.replace('-','')+"\n")
            infile.flush()
            cmd = "blastn -num_threads %i -task blastn -query '%s' -db '%s' -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send' -max_target_seqs %i" %\
                (num_threads, infile.name, db.sequences_fasta_file, max_target_seqs)
            logging.debug("Running cmd %s" % cmd)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            results_to_gather = []
            hit_counts = {}
            for line in iter(proc.stdout.readline,''):
                res = BlastQueryResultLine(line)
                query = queries[int(res.qseqid)]

                # Check we haven't hit max_target_seqs
                try:
                    hit_counts[res.qseqid] += 1
                except KeyError:
                    hit_counts[res.qseqid] = 1
                if hit_counts[res.qseqid] >= max_target_seqs:
                    logging.warn("The maximum number of target sequences returned by BLAST has been reached. Consider rerunning SingleM with an increased --max_hits cutoff.")

                query_length_original = len(query.sequence)
                query_length = len(query.sequence.replace('-',''))
                max_start = max([int(res.qstart),int(res.sstart)])-1
                pre_divergence = int(res.mismatch) + max_start

                # At this point, we do not know the length of the
                # subject sequence so we use only the query sequence
                # length, since the final divergence can only increase
                # when considering the subject sequence length.
                qtail_divergence = query_length-int(res.qend)
                divergence1 = pre_divergence + qtail_divergence
                logging.debug("Query %s hit of divergence1 %i" % (
                    res.qseqid, divergence1))
                if divergence1 <= max_divergence:
                    res.query = query
                    res.pre_divergence = pre_divergence
                    res.qtail_divergence = qtail_divergence
                    results_to_gather.append(res)

            logging.debug("Extracting %i sequences with blastdbcmd" % len(results_to_gather))
            # Extract all sequences in batch, to avoid repeated blastdbcmd calls.
            # Only extract one each result once
            dbseqs = db.extract_sequences_by_blast_ids(set([res.sseqid for res in results_to_gather]))
            sseqid_to_hits = {}
            for dbs in dbseqs:
                try:
                    sseqid_to_hits[dbs.sequence_id].add(dbs)
                except KeyError:
                    # It should be a set because two queries can hit the same subject
                    sseqid_to_hits[dbs.sequence_id] = set([dbs])

            # Calculate the final divergences with the extracted sequences.
            queries_subjects_divergences = []
            for res in results_to_gather:
                # All hits have the same sequence, so same divergence.
                subject_eg = next(iter(sseqid_to_hits[int(res.sseqid)]))
                subject_sequence = subject_eg.sequence
                # Simply align the sequences to avoid corner cases
                query_sequence = res.query.sequence
                if len(subject_sequence) != len(query_sequence):
                    raise Exception("At least for the moment, querying can only be carried out with 60bp OTU sequences, including gap characters.")
                divergence = 0
                for i, query_char in enumerate(query_sequence):
                    subject_char = subject_sequence[i]
                    if query_char != subject_char:
                        divergence = divergence + 1
                logging.debug("Query %s with subject '%s' had divergence2 %i" % (res.query.name, subject_eg.sequence_id, divergence))
                if divergence <= max_divergence:
                    queries_subjects_divergences.append([
                        res.query,
                        int(res.sseqid),
                        divergence])
            return BlastQueryResults(queries_subjects_divergences, sseqid_to_hits)

    def query_by_sqlite(self, queries, db_path):
        logging.info("Connecting to %s" % db_path)
        db = DatabaseManager({
        'sqlite3': {
            'driver': 'sqlite',
            'database': db_path
        }})
        Model.set_connection_resolver(db)

        results = []
        for query in queries:
            for entry in db.table('otus').where('sequence',seq).get():
                otu = OtuEntry()
                otu.marker = entry.marker
                otu.sample_name = entry.sample_name
                otu.sequence = entry.sequence
                otu.count = entry.count
                otu.coverage = entry.coverage
                otu.taxonomy = entry.taxonomy
                results.append(query, otu, 0)
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
