import subprocess
import string
from StringIO import StringIO

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
        cmd = "fxtract -XH -f /dev/stdin '%s'" % database_fasta_file

        process = subprocess.Popen(["bash", "-c", cmd],
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE)
        output, error = process.communicate('\n'.join(reads_to_extract))

        if process.returncode != 0:
            raise Exception("Extraction command '%s' failed with exitstatus %i" %
                            (cmd, process.returncode))

        seqs = []
        for name, seq, _ in SequenceIO().each(StringIO(output)):
            seqs.append(Sequence(name, seq))
        return seqs
