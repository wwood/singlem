import subprocess
import string
from io import StringIO
import extern

from graftm.sequence_io import Sequence, SequenceIO

class SequenceExtractor:
    '''Similar to graftm.SequenceExtractor except extracted sequences are read into
Python-land rather than output as a file.'''

    def extract_and_read(self, reads_to_extract, database_fasta_file):
        '''Extract the reads_to_extract from the database_fasta_file and return them.

        Parameters
        ----------
        reads_to_extract: Iterable of str
            IDs of reads to be extracted
        database_fasta_file: str
            path the fasta file that containing the reads

        Returns
        -------
        An array of graftm.sequence_io.Sequence objects'''
        cmd = "mfqe --output-uncompressed --fasta-read-name-lists /dev/stdin --input-fasta '{}' --output-fasta-files /dev/stdout".format(
            database_fasta_file)

        # Retrieve each sequence exactly once so mfqe does not croak
        output = extern.run(cmd, stdin='\n'.join(set(reads_to_extract)))

        seqs = []
        for name, seq, _ in SequenceIO().each(StringIO(output)):
            seqs.append(Sequence(name, seq))
        return seqs
