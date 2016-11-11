import os
import StringIO
from Bio import SeqIO

from otu_table import OtuTableEntry
import extern

class DBSequence(OtuTableEntry):
    sequence_id = None

    def fasta_defline(self):
        '''return the name of the sequence with all the info encoded'''
        return 'lcl|'+str(self.sequence_id)+" "+"|".join([
                                            self.marker,
                                            self.sample_name,
                                            str(self.count),
                                            self.taxonomy])

    def parse_from_fasta_define(self, defline):
        '''The opposite of fasta_define(), parse info into instance variables'''
        splits = defline.split(' ')
        if (splits) < 2: raise Exception("Parse exception: %s" % defline)
        self.sequence_id = int(splits[0].replace('lcl|',''))
        splits = ' '.join(splits[1:]).split('|')
        if len(splits) != 4: raise Exception("Parse exception 2: %s" % defline)

        self.marker, \
        self.sample_name,\
        self.count,\
        self.taxonomy = splits

        self.count = int(self.count)

class SequenceDatabase:
    version = 1

    @staticmethod
    def acquire(path):
        db = SequenceDatabase()
        db.sequences_fasta_file = os.path.join(path, "sequences.fasta")
        return db

    @staticmethod
    def create_from_otu_table(db_path, otu_table_collection):
        # ensure db does not already exist
        if os.path.exists(db_path):
            raise Exception("Cowardly refusing to overwrite already-existing database file '%s'" % db_path)
        os.makedirs(db_path)

        sequences_fasta_file = os.path.join(db_path, "sequences.fasta")
        with open(sequences_fasta_file, 'w') as fasta:
            sequence_id = 1
            for entry in otu_table_collection:
                dbseq = DBSequence()

                dbseq.marker = entry.marker
                dbseq.sample_name = entry.sample_name
                dbseq.sequence = entry.sequence
                dbseq.count = entry.count
                dbseq.taxonomy = entry.taxonomy
                dbseq.sequence_id = sequence_id

                fasta.write(">")
                fasta.write(dbseq.fasta_defline())
                fasta.write("\n")
                fasta.write(dbseq.sequence+"\n")

                sequence_id += 1

        cmd = "makeblastdb -in '%s' -dbtype nucl -parse_seqids" % sequences_fasta_file
        extern.run(cmd)

    def extract_sequence(self, sequence_id):
        cmd = "blastdbcmd -db '%s' -entry 'lcl|%s'" %\
            (self.sequences_fasta_file, sequence_id)
        stdout = extern.run(cmd)
        dbseq = DBSequence()
        for s in SeqIO.parse(StringIO.StringIO(stdout), "fasta"):
            if dbseq.sequence: raise Exception("Extracted multiple hits from sequence database, for sequence id %s" % sequence_id)
            dbseq.parse_from_fasta_define(s.description)
            dbseq.sequence = str(s.seq)
        return dbseq

    def extract_sequence_by_sseqid(self, sseqid):
        return self.extract_sequence(sseqid.replace('lcl|',''))
