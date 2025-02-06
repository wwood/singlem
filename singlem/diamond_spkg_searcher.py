import os
import logging
import extern

from subprocess import Popen, PIPE

from .utils import FastaNameToSampleName
from .run_via_os_system import run_via_os_system

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
            # running on reverse reads
            prefilter_dir = os.path.join(self._working_directory, 'prefilter_reverse')
        else:
            prefilter_dir = os.path.join(self._working_directory, 'prefilter_forward')
        os.mkdir(prefilter_dir)

        is_long_read = True
        
        for file in read_files:

            fasta_path = os.path.join(prefilter_dir,
                                      os.path.basename(file))
            if fasta_path[-3:] == '.gz':
                fasta_path = fasta_path[:-3] # remove .gz for destination files
            fasta_path = os.path.splitext(fasta_path)[0]+'.fna'
            
            f = open(fasta_path, 'w+') # create tempfile in working directory
            f.close()

            # #TODO: might get reid of the sed command and just use awk to print the first 3 columns
            # # this used to be under --db %s:  "| tee >(sed 's/^/>/; s/\\t/\\n/; s/\\t.*//' > %s) " \
            # cmd = "diamond blastx " \
            #       "--outfmt 6 qseqid full_qseq sseqid " \
            #       "--max-target-seqs 0 " \
            #       "--evalue 0.01 " \
            #       "%s " \
            #       "--threads %i " \
            #       "--query %s " \
            #       "--db %s " \
            #       "| awk '{print $1,$3,$2}'" % (
            #           performance_parameters,
            #           self._num_threads,
            #           file,
            #           diamond_database,
            #         #   fasta_path # this is the output file for the tee command
            #           )
            # logging.debug(cmd)

            # # Originially, we ran here via os.system rather than normal extern
            # # so reads can be piped in to singlem. However, this meant that
            # # errors and failed commands were ignored, sometimes causing
            # # successful return of singlem but an empty OTU table.
            # qseqid_sseqid = extern.run(cmd)
            # best_hits = {}

            # with open(fasta_path, 'a') as f:
            #     for line in qseqid_sseqid.splitlines():
            #         try:
            #             (qseqid,sseqid,full_qseq) = line.split(' ')

            #         except ValueError:
            #             raise Exception("Unexpected line format for DIAMOND output line '{}'".format(line))
                    
            #         if is_long_read:
            #             qseqid = qseqid + '~' + sseqid.split('~')[0]

            #         # TODO: idk how important that last check is
            #         # could maybe put in a bitscore check instead to ensure DIAMOND is outputing the best hits first
            #         if qseqid in best_hits and best_hits[qseqid] != sseqid:
            #             continue
            #         # else:
            #         #     if qseqid in best_hits and best_hits[qseqid] != sseqid:
            #         #         raise Exception("Multiple DIAMOND best hits detected for '{}'. This likely indicates that the input reads have non-unique names, possibly due to the same read appearing twice in a single input file".format(qseqid))

            #         best_hits[qseqid] = sseqid


            #         f.write('>' + qseqid + '\n')
            #         f.write(full_qseq + '\n')


            # cmd = [
            #     "diamond", "blastx",
            #     "--outfmt", "6", "qseqid", "full_qseq", "sseqid",
            #     "--max-target-seqs", "0",
            #     "--evalue", "0.01",
            #     "--threads", str(self._num_threads),
            #     "--query", file,
            #     "--db", diamond_database
            # ]

            # with range culling, etc
            cmd = [ 
                "diamond", "blastx",
                "--outfmt", "6", "qseqid", "full_qseq", "sseqid",
                "--max-target-seqs", "0",
                "--evalue", "0.01",
                "--frameshift", "15",
                "--range-culling",
                "-k", "1",
                "--threads", str(self._num_threads),
                "--query", file,
                "--db", diamond_database
            ]

            cmd.extend(performance_parameters.split())
            logging.debug(' '.join(cmd))

            best_hits = {}
            with Popen(cmd, stdout=PIPE, stderr=PIPE, text=True) as proc:
                with open(fasta_path, 'a') as fasta_file:
                    for line in proc.stdout:
                        try:
                            qseqid, full_qseq, sseqid = line.strip().split('\t')
                        except ValueError:
                            raise Exception(f"Unexpected line format for DIAMOND output line '{line.strip()}'")
                        
                        # TODO: potentially remove after testing, might not be necessary to explicitly check for long reads
                        if is_long_read:
                            # for forcing the simulated sshort reads to be treated as long reads
                            # if '-' in qseqid:
                            #     qseqid = qseqid.split('-')[0]

                            qseqid = qseqid + '~' + sseqid.split('~')[0]

                        if qseqid in best_hits:
                            continue

                        best_hits[qseqid] = sseqid

                        fasta_file.write(f'>{qseqid}\n{full_qseq}\n')    

                stderr_output = proc.stderr.read()
                if stderr_output:
                    logging.error(f"DIAMOND stderr: {stderr_output}")
                    raise Exception("DIAMOND failed")
            

            diamond_results.append(DiamondSearchResult(fasta_path, best_hits))

        return diamond_results

class DiamondSearchResult:
    def __init__(self, query_sequence_file, best_hits):
        self.query_sequences_file = query_sequence_file
        self.best_hits = best_hits

    def sample_name(self):
        return FastaNameToSampleName().fasta_to_name(self.query_sequences_file)