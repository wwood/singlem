import subprocess

class SequenceExtractor:
    def extract(self, reads_to_extract, database_fasta_file, output_file):
        '''Extract the reads_to_extract from the database_fasta_file and put them in
        output_file.

        Parameters
        ----------
        reads_to_extract: Iterable of str
            IDs of reads to be extracted
        database_fasta_file: str
            path the fasta file that containing the reads
        output_file: str
            path to the file where they are put
        
        Returns
        -------
        Nothing'''
        cmd = "fxtract -XH -f /dev/stdin '%s' > %s" % (
            database_fasta_file, output_file)

        process = subprocess.Popen(["bash", "-c", cmd], 
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE)
        process.communicate('\n'.join(reads_to_extract))

    def extract_and_reverse_complement(self, reads_to_extract,
                                       database_fasta_file, output_file):
        '''As per extract except also reverse complement the sequences.'''
        cmd = "fxtract -XH -f /dev/stdin '%s' | " % database_fasta_file + \
              "seqmagick convert --reverse-complement --input-format fasta "\
              "--output-format fasta - %s" % (output_file)

        process = subprocess.Popen(["bash", "-c", cmd], 
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE)
        process.communicate('\n'.join(reads_to_extract))
        
