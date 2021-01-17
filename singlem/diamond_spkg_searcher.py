import os
import extern

class DiamondSpkgSearcher:
    def __init__(self, num_threads, working_directory):
        self._num_threads = num_threads
        self._working_directory = working_directory

    def run_diamond(self, hmms, forward_read_files, reverse_read_files):
        for pkg in hmms:
            if not pkg.is_protein_package():
                raise Exception(
                    "DIAMOND prefilter cannot be used with nucleotide SingleM packages")
        
        dmnd = hmms.get_dmnd()
        new_forward_read_files = self._prefilter(dmnd, forward_read_files, False)
        new_reverse_read_files = None
        if reverse_read_files != None:
            new_reverse_read_files = self._prefilter(dmnd, reverse_read_files, True)

        return (new_forward_read_files, new_reverse_read_files)

    def _prefilter(self, diamond_database, read_files, is_reverse_reads):
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
        path to fasta file of filtered reads
        '''
        filtered_reads = []
        if is_reverse_reads:
            # running on reverse reads
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
            
            cmd = "zcat -f %s | " \
                  "diamond blastx " \
                  "--outfmt 6 qseqid full_qseq " \
                  "--max-target-seqs 1 " \
                  "--evalue 0.01 " \
                  "--index-chunks 1 " \
                  "--threads %i " \
                  "--query - " \
                  "--db %s " \
                  "| sed -e 's/^/>/' -e 's/\\t/\\n/' > %s" % (
                      file,
                      self._num_threads,
                      diamond_database,
                      fasta_path)
            extern.run(cmd)
            filtered_reads.append(fasta_path)
            
        return filtered_reads