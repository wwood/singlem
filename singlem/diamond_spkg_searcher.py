import os
import logging

from subprocess import Popen, PIPE

from .utils import FastaNameToSampleName

class DiamondSpkgSearcher:
    def __init__(self, num_threads, working_directory, context_window=None):
        self._num_threads = num_threads
        self._working_directory = working_directory
        self._context_window = context_window

    def run_diamond(self, hmms, forward_read_files, reverse_read_files, performance_parameters, diamond_db, min_orf_length):
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

        # Calculate sample names here so they can be used in the reverse read
        # processing too. If this does not happen, getting the full length
        # sequences is hard later on because the reverse read basename is not
        # necessarily the sample as the forward read one.
        sample_names = []
        for file in forward_read_files:
            basename = os.path.basename(file)
            if basename[-3:] == '.gz':
                basename = basename[:-3] # remove .gz for destination files
            basename = os.path.splitext(basename)[0]+'.fna'
            sample_names.append(basename)

        fwds = self._prefilter(dmnd, forward_read_files, False, performance_parameters, sample_names, min_orf_length)
        revs = None
        if reverse_read_files != None:
            revs = self._prefilter(dmnd, reverse_read_files, True, performance_parameters, sample_names, min_orf_length)

        return (fwds, revs)

    def _prefilter(self, diamond_database, read_files, is_reverse_reads, performance_parameters, sample_names=None, min_orf_length=None):
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

        if sample_names is None:
            sample_names = []
            for file in read_files:
                basename = os.path.basename(file)
                if basename[-3:] == '.gz':
                    basename = basename[:-3]
                sample_names.append(os.path.splitext(basename)[0]+'.fna')
        if min_orf_length is None:
            min_orf_length = 0

        if is_reverse_reads:
            prefilter_dir = os.path.join(self._working_directory, 'prefilter_reverse')
        else:
            prefilter_dir = os.path.join(self._working_directory, 'prefilter_forward')
        os.mkdir(prefilter_dir)

        for (file, sample_name) in zip(read_files, sample_names):
            fasta_path = os.path.join(prefilter_dir, sample_name)
            
            # TODO: Why is this w+ needed? I suppose it does no harm though so
            # leaving it for the moment.
            f = open(fasta_path, 'w+') # create tempfile in working directory
            f.close()
            previous_sample_name = os.path.basename(fasta_path)
            full_qseq_fasta_path = fasta_path + '.full_qseqs'
            full_qseq_f = open(full_qseq_fasta_path, 'w')
            full_qseq_f.close()

            # DIAMOND command now with range culling, etc
            #
            # Note these parameters should align with those used in the taxonomy
            # search step later on in pipe, otherwise some reads that pass the
            # prefilter will not be assigned any taxonomy.
            cmd = [ 
                "diamond", "blastx",
                "--outfmt", "6", "qseqid", "full_qseq", "sseqid", "qstart", "qend",
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
            query_sequence_lengths = {}
            alignment_positions = {}
            # using Popen to stream the output
            with Popen(cmd, stdout=PIPE, stderr=PIPE, text=True) as proc:
                seen_full_qseqs = set()
                with open(fasta_path, 'a') as fasta_file, open(full_qseq_fasta_path, 'a') as full_qseq_f:
                    for line in proc.stdout:
                        try:
                            parts = line.strip().split('\t')
                            if len(parts) == 5:
                                qseqid, full_qseq, sseqid, qstart, qend = parts
                                query_length = len(full_qseq)
                            elif len(parts) == 6:
                                qseqid, full_qseq, sseqid, qstart, qend, qlen = parts
                                query_length = int(qlen)
                            else:
                                raise ValueError
                        except ValueError:
                            raise Exception(f"Unexpected line format for DIAMOND output line '{line.strip()}'")
                    
                        aligned_start_raw = min(int(qstart), int(qend))
                        aligned_end_raw = max(int(qstart), int(qend))

                        # creating new read index to account for multiple hits
                        # by concating the read_name with the marker_gene_name, we can ensure only 1 gene copy per read
                        # TODO: add an option to let all unique genes through with range-uclling 
                        if self._context_window is not None:
                            context_start = max(1, aligned_start_raw - self._context_window)
                            context_end = min(query_length, aligned_end_raw + self._context_window)
                            unique_qseqid = f"{qseqid}:{context_start}-{context_end},{query_length}••{sseqid.split('~')[0]}"
                            alignment_start = context_start
                            alignment_end = context_end
                        else:
                            unique_qseqid = qseqid + '••' + sseqid.split('~')[0]
                            alignment_start = aligned_start_raw
                            alignment_end = aligned_end_raw

                        # extra check to make sure we're not overwriting a better hit
                        if unique_qseqid in best_hits:
                            continue

                        # store the best hit and sequence length for each query sequence to feed into the next steps
                        best_hits[unique_qseqid] = sseqid
                        query_sequence_lengths[unique_qseqid] = query_length
                        if qseqid not in alignment_positions:
                            alignment_positions[qseqid] = (alignment_start, alignment_end, query_length)

                        # Only write the part of the query sequence that aligned
                        # to the prefilter database to the fasta file.
                        # qstat can be > qend if the sequences is reverse complemented.
                        if self._context_window is not None:
                            qstart_int = alignment_start
                            qend_int = alignment_end
                        else:
                            qstart_int = alignment_start
                            qend_int = alignment_end
                            # To ensure that that the ORF finder later on finds something >=72bp, extend the range
                            # either side if possible.
                            extra_bp = min_orf_length
                            qstart_int = max(1, qstart_int - extra_bp)
                            qend_int = min(query_length, qend_int + extra_bp)
                        logging.debug(f"DIAMOND hit: {qseqid} {sseqid} {qstart} {qend} => {qstart_int} {qend_int}")
                        truncated_qseq = full_qseq[int(qstart_int)-1:int(qend_int)]
                        fasta_file.write(f'>{unique_qseqid}\n{truncated_qseq}\n')

                        # Write the full read sequence to the fasta file as well
                        # so it can be included in the archive OTU table later
                        # on. Write without the unique_qseqid modification, to save space.
                        if qseqid not in seen_full_qseqs:
                            if self._context_window is not None:
                                full_qseq_f.write(f'>{qseqid}\n{truncated_qseq}\n')
                            else:
                                full_qseq_f.write(f'>{qseqid}\n{full_qseq}\n')
                            seen_full_qseqs.add(qseqid)


                # check for DIAMOND errors
                stderr_output = proc.stderr.read()
                if stderr_output:
                    logging.error(f"DIAMOND stderr: {stderr_output}")
                    raise Exception("DIAMOND failed")
                
                # check for non-zero return code
                return_code = proc.wait()
                if return_code != 0:
                    raise Exception(f"DIAMOND failed with return code {return_code}, but no stderr output")

            diamond_results.append(DiamondSearchResult(
                fasta_path,
                full_qseq_fasta_path,
                best_hits,
                query_sequence_lengths,
                alignment_positions))

        return diamond_results

class DiamondSearchResult:
    def __init__(self, query_sequence_file, full_query_sequences_file, best_hits, query_sequence_lengths, alignment_positions):
        self.query_sequences_file = query_sequence_file
        self.full_query_sequences_file = full_query_sequences_file
        self.best_hits = best_hits
        self.query_sequence_lengths = query_sequence_lengths
        self.alignment_positions = alignment_positions

    def sample_name(self):
        return FastaNameToSampleName().fasta_to_name(self.query_sequences_file)
