import tempfile
import extern
import logging
import subprocess
import sys

from sequence_database import SequenceDatabase
from sequence_classes import SeqReader
from query_formatters import NameSequenceQueryDefinition, NamedQueryDefinition
from query_formatters import SparseResultFormatter, DenseResultFormatter
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
            query_names = ['unnamed_sequence']
            query_sequences = [query_sequence]
        elif query_otu_table:
            query_names = []
            query_sequences = []
            otus = OtuTableCollection()
            otus.add_otu_table(open(query_otu_table))
            for e in otus:
                query_sequences.append(e.sequence)
                query_names.append(';'.join([e.sample_name,e.marker]))
        elif query_fasta:
            query_names = []
            query_sequences = []
            for name, seq, _ in SeqReader().readfq(open(query_fasta)):
                query_names.append(name)
                query_sequences.append(seq)
        else:
            raise Exception("No query option specified, cannot continue")

        # blast the query against the database, output as csv
        found_distances_and_names = []
        with tempfile.NamedTemporaryFile(prefix='singlem_query') as infile:
            for i, sequence in enumerate(query_sequences):
                infile.write(">%i\n" % i)
                infile.write(sequence.replace('-','')+"\n")
            infile.flush()

            cmd = "blastn -num_threads %i -task blastn -query '%s' -db '%s' -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send' -max_target_seqs %i" %\
                (num_threads, infile.name, db.sequences_fasta_file, max_target_seqs)
            logging.debug("Running cmd %s" % cmd)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            last_index = None
            last_differences_and_names = []
            results_to_gather = []
            for line in iter(proc.stdout.readline,''):
                res = QueryResultLine(line)
                index = int(res.qseqid)
                #TODO: check we haven't come up against the max_target_seqs barrier
                query_length_original = len(query_sequences[index])
                query_length = len(query_sequences[index].replace('-',''))
                max_start = max([int(res.qstart),int(res.sstart)])-1
                pre_divergence = int(res.mismatch) + max_start
                # At this point, we do not know the length of the subject sequence so we use only the query sequence length, since the final divergence can only increase when considering the subject sequence length.
                qtail_divergence = query_length-int(res.qend)
                divergence1 = pre_divergence + qtail_divergence
                if divergence1 <= max_divergence:
                    res.query_index = index
                    res.pre_divergence = pre_divergence
                    res.qtail_divergence = qtail_divergence
                    results_to_gather.append(res)

            # Extract all sequences in batch, to avoid repeated blastdbcmd calls.
            dbseqs = db.extract_sequences_by_blast_ids([res.sseqid for res in results_to_gather])
            # Calculate the final divergences with the extracted sequences.
            found_query_names = []
            found_query_sequences = []
            last_query_index = None
            divergences_and_subjects = []
            for i, res in enumerate(results_to_gather):
                subject = dbseqs[i]
                # Simply align the sequences to avoid corner cases
                query_sequence = query_sequences[res.query_index]
                if len(subject.sequence) != len(query_sequence):
                    raise Exception("At least for the moment, querying can only be carried out with 60bp OTU sequences, including gap characters.")
                divergence = 0
                for i, query_char in enumerate(query_sequence):
                    if query_char != subject.sequence[i]:
                        divergence = divergence + 1
                if divergence <= max_divergence:
                    if res.query_index != last_query_index:
                        last_query_index = res.query_index
                        found_query_names.append(query_names[res.query_index])
                        found_query_sequences.append(query_sequences[res.query_index])
                        divergences_and_subjects.append([[divergence, subject]])
                    else:
                        divergences_and_subjects[-1].append([divergence, subject])

        if query_fasta:
            namedef = NamedQueryDefinition(found_query_names)
        else:
            namedef = NameSequenceQueryDefinition(found_query_names, found_query_sequences)

        if output_style == 'sparse':
            formatter = SparseResultFormatter(namedef, divergences_and_subjects)
        elif output_style == 'dense':
            formatter = DenseResultFormatter(namedef, divergences_and_subjects)
        else:
            raise Exception()
        formatter.write(sys.stdout)


class QueryResultLine:
    def __init__(self, blast_output_line):
        self.qseqid, self.sseqid, _, _, self.mismatch, self.gaps, self.qstart,\
            self.qend, self.sstart, \
            self.send = blast_output_line.strip().split("\t")[:10]
