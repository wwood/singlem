import os
import logging

from subprocess import Popen, PIPE

from .utils import FastaNameToSampleName

class DiamondSpkgSearcher:
    def __init__(self, num_threads, working_directory):
        self._num_threads = num_threads
        self._working_directory = working_directory

    def run_diamond(self, hmms, forward_read_files, reverse_read_files, performance_parameters, diamond_db):
        '''Run a single DIAMOND run for each of the forward_read_files against a 
        combined database of all sequences from the singlem package set given.

        diamond_db: None or str
            path to provide to DIAMOND for its DB (or None to create on the fly)

        Returns
        -------
        (fwds, revs) where fwds is a list of DiamondSearchResult objects, and revs is
        the same, or None if reverse read input is None
        '''
        for pkg in hmms:
            if not pkg.is_protein_package():
                raise Exception(
                    "DIAMOND prefilter cannot be used with nucleotide SingleM packages")
        
        if diamond_db is None:
            dmnd = hmms.get_dmnd()
        else:
            dmnd = diamond_db
        fwds = self._prefilter(dmnd, forward_read_files, False, performance_parameters)
        revs = None
        if reverse_read_files != None:
            revs = self._prefilter(dmnd, reverse_read_files, True, performance_parameters)

        return (fwds, revs)

    def _prefilter(self, diamond_database, read_files, is_reverse_reads, performance_parameters):
        '''Find all reads that match the DIAMOND database in the 
        singlem_package database.
        Parameters
        ----------
        diamond_database: dmnd 
            DIAMOND database from the SingleM packages to search reads with
        read_files: list of str 
            paths to the sequences to be searched
        reverse: boolean
            check if using reverse reads
        Returns
        -------
        Array of DiamondSearchResult objects, one for each read file
        '''
        diamond_results = []

        if is_reverse_reads:
            prefilter_dir = os.path.join(self._working_directory, 'prefilter_reverse')
        else:
            prefilter_dir = os.path.join(self._working_directory, 'prefilter_forward')
        os.mkdir(prefilter_dir)

        
        for file in read_files:

            fasta_path = os.path.join(prefilter_dir,
                                      os.path.basename(file))
            if fasta_path[-3:] == '.gz':
                fasta_path = fasta_path[:-3] # remove .gz for destination files
            fasta_path = os.path.splitext(fasta_path)[0]+'.fna'
            
            f = open(fasta_path, 'w+') # create tempfile in working directory
            f.close()

            # DIAMOND command
            # now with range culling, etc
            cmd = [ 
                "diamond", "blastx",
                "--outfmt", "6", "qseqid", "full_qseq", "sseqid", "qstart",
                "--max-target-seqs", "1",
                "--evalue", "0.01",
                "--frameshift", "15",
                "--range-culling",
                "--range-cover", "1",
                "--threads", str(self._num_threads),
                "--query", file,
                "--db", diamond_database
            ]

            cmd.extend(performance_parameters.split())
            logging.debug(' '.join(cmd))

            best_hits = {}
            # using Popen to stream the output
            with Popen(cmd, stdout=PIPE, stderr=PIPE, text=True) as proc:
                with open(fasta_path, 'a') as fasta_file:
                    for line in proc.stdout:
                        try:
                            qseqid, full_qseq, sseqid, qstart = line.strip().split('\t')
                        except ValueError:
                            raise Exception(f"Unexpected line format for DIAMOND output line '{line.strip()}'")
                    
                        # creating new read index to account for multiple hits
                        # by concating the read_name with the marker_gene_name, we can ensure only 1 gene copy per read
                        # TODO: add an option to let all unique genes through with range-uclling 
                        qseqid = qseqid + '••' + sseqid.split('~')[0]

                        # extra check to make sure we're not overwriting a better hit
                        if qseqid in best_hits:
                            continue

                        # store the best hit for each query sequence to feed into the next steps
                        best_hits[qseqid] = sseqid

                        # write the query sequence to a file
                        fasta_file.write(f'>{qseqid}\n{full_qseq}\n')    

                # check for DIAMOND errors
                stderr_output = proc.stderr.read()
                if stderr_output:
                    logging.error(f"DIAMOND stderr: {stderr_output}")
                    raise Exception("DIAMOND failed")
                
                # check for non-zero return code
                return_code = proc.wait()
                if return_code != 0:
                    raise Exception(f"DIAMOND failed with return code {return_code}, but no stderr output")

            diamond_results.append(DiamondSearchResult(fasta_path, best_hits))

        return diamond_results

class DiamondSearchResult:
    def __init__(self, query_sequence_file, best_hits):
        self.query_sequences_file = query_sequence_file
        self.best_hits = best_hits

    def sample_name(self):
        return FastaNameToSampleName().fasta_to_name(self.query_sequences_file)