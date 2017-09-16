import csv
import os
import tempfile
import StringIO
import logging
import subprocess
import itertools
import sqlite3
import glob

from itertools import izip_longest
from threading import Thread
from Bio import SeqIO

from otu_table import OtuTableEntry
import extern

class DBSequence(OtuTableEntry):
    sequence_id = None

    def fasta_defline(self):
        '''return the name of the sequence with all the info encoded'''
        return "|".join([
            self.marker,
            self.sample_name,
            str(self.count),
            self.taxonomy])

    @staticmethod
    def parse_from_fasta_define(defline):
        '''The opposite of fasta_define(), parse info into instance variables in
        an array of DBSequence objects
        '''
        splits = defline.split(' ')
        if (splits) < 2: raise Exception("Parse exception: %s" % defline)
        sequence_id = int(splits[0].replace('lcl|',''))

        splits = ' '.join(splits[1:])
        to_return = []
        for subsplit in splits.split(SequenceDatabase.DEFLINE_DELIMITER_CHARACTER):
            splits2 = subsplit.split('|')
            if len(splits2) != 4: raise Exception("Parse exception 2: %s" % defline)
            dbseq = DBSequence()
            dbseq.sequence_id = sequence_id
            dbseq.marker,\
            dbseq.sample_name,\
            dbseq.count,\
            dbseq.taxonomy = splits2
            dbseq.count = int(dbseq.count)
            to_return.append(dbseq)
        return to_return

class SequenceDatabase:
    version = 2
    GAP_REPLACEMENT_CHARACTER = 'Y'
    DEFLINE_DELIMITER_CHARACTER = '~'
    SQLITE_DB_NAME = 'otus.sqlite3'
    _marker_to_smafadb = {}

    def add_smafa_db(self, marker_name, smafa_db):
        self._marker_to_smafadb[marker_name] = smafa_db

    def get_smafa_db(self, marker_name):
        if marker_name in self._marker_to_smafadb:
            return self._marker_to_smafadb[marker_name]
        else:
            logging.debug("No smafa DB found for %s" % marker_name)
            return None

    def smafa_dbs(self):
        return self._marker_to_smafadb.values()

    @staticmethod
    def acquire(path):
        db = SequenceDatabase()
        db.sqlite_file = os.path.join(path, SequenceDatabase.SQLITE_DB_NAME)
        smafas = glob.glob("%s/*.smafadb" % path)
        logging.debug("Found smafadbs: %s" % ", ".join(smafas))
        if len(smafas) == 0: raise Exception("No smafa DBs found in DB")
        for g in smafas:
            marker = os.path.basename(g).replace('.smafadb','') # TODO: Make a contents file.
            db.add_smafa_db(marker, g)
        return db

    @staticmethod
    def grouper(iterable, n):
        args = [iter(iterable)] * n
        return izip_longest(*args, fillvalue=None)

    @staticmethod
    def create_from_otu_table(db_path, otu_table_collection):
        # ensure db does not already exist
        if os.path.exists(db_path):
            raise Exception("Cowardly refusing to overwrite already-existing database file '%s'" % db_path)
        os.makedirs(db_path)

        # setup sqlite DB
        sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
        logging.debug("Connecting to db %s" % sqlite_db_path)
        db = sqlite3.connect(sqlite_db_path)
        c = db.cursor()
        c.execute("CREATE TABLE otus (marker text, sample_name text,"
                  " sequence text, num_hits int, coverage float, taxonomy text)")
        db.commit()

        gene_to_tempfile = {}

        # logging.info("Writing FASTA file to %s .." % sequences_fasta_file)
        # cmd = "sort -S 20%" # The default sort buffer is too low IMO, use a larger one.
        # sorter = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE, stdout=subprocess.PIPE)
        # fasta_writing_thread = Thread(target=SequenceDatabase.write_dereplicated_fasta_file,
        #                               args=[sorter.stdout, fasta])
        # fasta_writing_thread.start()
        chunksize = 10000 # Run in chunks for sqlite insert performance.
        for chunk in SequenceDatabase.grouper(otu_table_collection, chunksize):
            chunk_list = []
            for entry in chunk:
                if entry is not None: # Is None when padded in last chunk.
                    chunk_list.append((entry.marker, entry.sample_name, entry.sequence, entry.count,
                        entry.coverage, entry.taxonomy))
                    dbseq = DBSequence()
                    dbseq.marker = entry.marker
                    dbseq.sample_name = entry.sample_name
                    dbseq.sequence = entry.sequence
                    dbseq.count = entry.count
                    dbseq.taxonomy = entry.taxonomy

                    if entry.marker not in gene_to_tempfile:
                        gene_to_tempfile[entry.marker] = tempfile.NamedTemporaryFile(prefix='singlem-makedb')
                    tf = gene_to_tempfile[entry.marker]
                    tf.write("%s\n" % entry.sequence)

            c.executemany("INSERT INTO otus(marker, sample_name, sequence, num_hits, "
                          "coverage, taxonomy) VALUES(?,?,?,?,?,?)",
                          chunk_list)

        logging.info("Creating SQLite indices")
        c.execute("CREATE INDEX otu_sequence on otus (sequence)")
        c.execute("CREATE INDEX otu_sample_name on otus (sample_name)")
        db.commit()

        # Run smafa on each of the genomes
        logging.info("Running smafas")
        for marker_name, tf in gene_to_tempfile.items():
            smafa = "%s.smafadb" % os.path.join(db_path, marker_name)
            tf.flush()
            cmd = "sort -S20%% '%s' |uniq |awk '{print \">1\\n\" $1}' |smafa makedb /dev/stdin '%s'" % (
                tf.name, smafa)
            logging.debug("Running cmd: %s", cmd)
            logging.info("Formatting smafa database %s .." % smafa)
            extern.run(cmd)
            tf.close()
        logging.info("Finished")
