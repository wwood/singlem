import csv
import os
import tempfile
import StringIO
import logging
import subprocess
import itertools
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
    version = 1
    GAP_REPLACEMENT_CHARACTER = 'Y'
    DEFLINE_DELIMITER_CHARACTER = '~'

    @staticmethod
    def acquire(path):
        db = SequenceDatabase()
        db.sequences_fasta_file = os.path.join(path, "sequences.fasta")
        return db

    @staticmethod
    def write_dereplicated_fasta_file(input_stream, output_stream):
        logging.debug("Starting to write FASTA file in thread")
        last_sequence = None
        ids_to_print = []
        num_printed = 0
        sequence_id = 1
        for d in csv.reader(input_stream, delimiter="\t"):
            if SequenceDatabase.DEFLINE_DELIMITER_CHARACTER in d[1]:
                raise Exception(
                    "Taxonomy and sample names cannot have '%s' in them" % \
                    SequenceDatabase.DEFLINE_DELIMITER_CHARACTER)
            if d[0] == last_sequence:
                ids_to_print.append(d[1])
            else:
                if ids_to_print != []: # if not the first line
                    output_stream.write(">lcl|%i " % sequence_id)
                    sequence_id += 1
                    output_stream.write(SequenceDatabase.DEFLINE_DELIMITER_CHARACTER.join(ids_to_print))
                    output_stream.write("\n")
                    output_stream.write(last_sequence)
                    output_stream.write("\n")
                ids_to_print = [d[1]]
            last_sequence = d[0]
        if ids_to_print != []:
            # Print the last sequence
            output_stream.write(">lcl|%i " % sequence_id)
            output_stream.write(SequenceDatabase.DEFLINE_DELIMITER_CHARACTER.join(ids_to_print))
            output_stream.write("\n")
            output_stream.write(d[0])
            output_stream.write("\n")
        logging.debug("Finished writing fasta in thread")

    @staticmethod
    def create_from_otu_table(db_path, otu_table_collection):
        # ensure db does not already exist
        if os.path.exists(db_path):
            raise Exception("Cowardly refusing to overwrite already-existing database file '%s'" % db_path)
        os.makedirs(db_path)

        # 100% cluster the database.
        sequences_fasta_file = os.path.join(db_path, "sequences.fasta")
        with open(sequences_fasta_file, 'w') as fasta:
            logging.info("Writing FASTA file to %s .." % sequences_fasta_file)
            cmd = "sort -S 20%" # The default sort buffer is too low IMO, use a larger one.
            sorter = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE, stdout=subprocess.PIPE)
            fasta_writing_thread = Thread(target=SequenceDatabase.write_dereplicated_fasta_file,
                                          args=[sorter.stdout, fasta])
            fasta_writing_thread.start()
            for entry in otu_table_collection:
                dbseq = DBSequence()
                dbseq.marker = entry.marker
                dbseq.sample_name = entry.sample_name
                dbseq.sequence = entry.sequence
                dbseq.count = entry.count
                dbseq.taxonomy = entry.taxonomy

                # Replace gaps with a replacement char so we can tell the difference between - and N.
                if SequenceDatabase.GAP_REPLACEMENT_CHARACTER in dbseq.sequence:
                    logging.warn("Attempting to create database with the reserved character %s, calculated divergences may be slightly incorrect" % SequenceDatabase.GAP_REPLACEMENT_CHARACTER)
                sorter.stdin.write(dbseq.sequence.replace('-',SequenceDatabase.GAP_REPLACEMENT_CHARACTER)+"\t")
                sorter.stdin.write(dbseq.fasta_defline() + "\n")

            sorter.stdin.close()
            fasta_writing_thread.join()

        cmd = "makeblastdb -in '%s' -dbtype nucl -parse_seqids" % sequences_fasta_file
        logging.info("Generating BLAST database ..")
        extern.run(cmd)
        logging.info("Finished")

    def extract_sequences(self, sequence_id):
        cmd = "blastdbcmd -db '%s' -entry 'lcl|%s'" %\
            (self.sequences_fasta_file, sequence_id)
        stdout = extern.run(cmd)
        for s in SeqIO.parse(StringIO.StringIO(stdout), "fasta"):
            if dbseq.sequence: raise Exception("Extracted multiple hits from sequence database, for sequence id %s" % sequence_id)
            dbseq.parse_from_fasta_define(s.description)
            dbseq.sequence = str(s.seq)
        return dbseq

    def extract_sequences_by_blast_ids(self, blast_ids):
        num_blast_ids = 0
        with tempfile.NamedTemporaryFile(prefix='singlem_query_for_blastdbcmd') as batchfile:
            for bid in blast_ids:
                batchfile.write(bid)
                batchfile.write("\n")
                num_blast_ids += 1
            batchfile.flush()
            cmd = "blastdbcmd -db '%s' -entry_batch '%s'"\
                   % (self.sequences_fasta_file, batchfile.name)
            stdout = extern.run(cmd)
            dbseqs = []
            for s in SeqIO.parse(StringIO.StringIO(stdout), "fasta"):
                new_seqs = DBSequence.parse_from_fasta_define(s.description)
                for single_dbseq in new_seqs:
                    single_dbseq.sequence = str(s.seq)
                dbseqs.append(new_seqs)
            if len(dbseqs) != num_blast_ids:
                raise Exception("Unexpected number of returned sequences from blastdbcmd")
        return list(itertools.chain(*dbseqs))
