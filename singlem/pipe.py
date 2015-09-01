import tempdir
import logging
import os
import shutil
import extern

from singlem import MetagenomeOtuFinder, \
    HmmDatabase, TaxonomyFile, SeqReader
import itertools
import tempfile
from known_otu_table import KnownOtuTable
from string import split
import subprocess
    
class SearchPipe:
    def run(self, **kwargs):
        forward_read_files = kwargs.pop('sequences')
        output_otu_table = kwargs.pop('otu_table')
        num_threads = kwargs.pop('threads')
        bootstrap_contigs = kwargs.pop('bootstrap_contigs')
        known_otu_tables = kwargs.pop('known_otu_tables')
        graftm_assignment_method = kwargs.pop('assignment_method')
        output_extras = kwargs.pop('output_extras')
        
        working_directory = kwargs.pop('working_directory')
        force = kwargs.pop('force')
        previous_graftm_search_directory = kwargs.pop('previous_graftm_search_directory')
        previous_graftm_separate_directory = kwargs.pop('previous_graftm_separate_directory')
        
        hmms = HmmDatabase()
        
        using_temporary_working_directory = working_directory is None
        if using_temporary_working_directory:
            tmp = tempdir.TempDir()
            working_directory = tmp.name
        else:
            working_directory = working_directory
            if os.path.exists(working_directory):
                if previous_graftm_search_directory:
                    logging.debug("Using previously half-run singlem results")
                elif force:
                    logging.info("Overwriting directory %s" % working_directory)
                    shutil.rmtree(working_directory)
                    os.mkdir(working_directory)
                else:
                    raise Exception("Working directory '%s' already exists, not continuing" % working_directory)
            else:
                os.mkdir(working_directory)
        logging.debug("Using working directory %s" % working_directory)
        
        def return_cleanly():
            if using_temporary_working_directory: tmp.dissolve()
            
        if bootstrap_contigs:
            # Create HMMs from each of the search HMMs based on the given contigs
            bootstrap_hmms = {}
            logging.info("Generating HMMs by bootstrap..")
            # TOOD: this can be combined into a single run over the data with
            # orfm, this is a bit wasteful going over it 3 times. Instead
            # best to interface with graftm directly, I guess.
            for hmm in hmms:
                bootstrap_hmm_file = os.path.join(working_directory,
                              "bootstrap_%s" % os.path.basename(hmm.hmm_filename))
                bootstrap_hmms[hmm.hmm_filename] = bootstrap_hmm_file
                logging.debug("Finding bootstrap hmm %s" % bootstrap_hmm_file)
                cmd = "graftM bootstrap --contigs %s --search_hmm_files %s "\
                    " --verbosity 2 --output_hmm %s" %(\
                      ' '.join(bootstrap_contigs),
                      hmm.hmm_path(),
                      bootstrap_hmm_file)
                logging.debug("Running cmd: %s" % cmd)
                extern.run(cmd)
        
        if not previous_graftm_search_directory:
            graftm_search_directory = os.path.join(working_directory, 'graftm_search')
        
            # run graftm across all the HMMs
            logging.info("Using as input %i different forward read sets e.g. %s" % (len(forward_read_files),
                                                                                forward_read_files[0]))
            cmd = "graftM graft --threads %i --forward %s "\
                "--search_hmm_files %s --search_and_align_only "\
                "--output_directory %s --aln_hmm_file %s --verbosity 2 "\
                "--input_sequence_type nucleotide"\
                                 % (num_threads,
                                    ' '.join(forward_read_files),
                                    ' '.join(hmms.hmm_paths()),
                                    graftm_search_directory,
                                    hmms.hmm_paths()[0])
            if bootstrap_contigs:
                cmd += " --search_hmm_files %s" % ' '.join(
                    itertools.chain(
                        (f for f in bootstrap_hmms.values() if os.path.isfile(f)),
                        hmms.hmm_paths()))
            logging.info("Running GraftM to find particular reads..")
            extern.run(cmd)
            logging.info("Finished running GraftM search phase")
        else:
            graftm_search_directory = previous_graftm_search_directory
            logging.info("Using existing graftM directory %s" % previous_graftm_search_directory)
    
        # Get the names of the samples from the graftm directory
        sample_names = [f for f in os.listdir(graftm_search_directory) \
                        if os.path.isdir(os.path.join(graftm_search_directory, f))]
        logging.debug("Recovered %i samples from graftm search output e.g. %s" \
                     % (len(sample_names), sample_names[0]))
        
        if previous_graftm_separate_directory:
            graftm_separate_directory_base = previous_graftm_separate_directory
        else:
            graftm_separate_directory_base = os.path.join(working_directory, 'graftm_separates')
            os.mkdir(graftm_separate_directory_base)
            
            with tempfile.NamedTemporaryFile() as samples_file:
                # Don't look at samples that have no hits
                viable_sample_names = []
                for sample_name in sample_names:
                    if os.stat("%s/%s/%s_hits.fa" % (graftm_search_directory,
                                                  sample_name,
                                                  sample_name)).st_size > 0:
                        viable_sample_names.append(sample_name)
                    else:
                        logging.warn("Sample '%s' does not appear to contain any hits" % sample_name)
                if len(viable_sample_names) == 0:
                    logging.info("No reads identified in any samples, stopping")
                    return_cleanly()
                    return
                else:
                    logging.debug("Found %i samples with reads identified" % 
                                  len(viable_sample_names))
                samples_file.write("\n".join(viable_sample_names))
                samples_file.flush()
                
                logging.info("Running separate alignments in GraftM..")
                commands = []
                for sample_name in viable_sample_names:
                    for hmm in hmms:
                        cmd = "graftM graft --threads %i --verbosity 2 "\
                             "--forward %s/%s/%s_hits.fa "\
                             "--graftm_package %s --output_directory %s/%s_vs_%s "\
                             "--input_sequence_type nucleotide "\
                             "--search_and_align_only" % (\
                                    1, #use 1 thread since most likely better to parallelise processes with extern, not threads here
                                    graftm_search_directory,
                                    sample_name,
                                    sample_name,
                                    hmm.gpkg_path,
                                    graftm_separate_directory_base,
                                    sample_name,
                                    os.path.basename(hmm.gpkg_path))
                        if bootstrap_contigs:
                            bootstrap_hmm = bootstrap_hmms[hmm.hmm_filename]
                            if os.path.isfile(bootstrap_hmm):
                                cmd += " --search_hmm_files %s %s" % (
                                            bootstrap_hmm,
                                            hmm.hmm_path())
                        commands.append(cmd)
                extern.run_many(commands, num_threads=num_threads)
                
        sample_to_gpkg_to_input_sequences = self._extract_relevant_reads(graftm_separate_directory_base, sample_names, hmms)
        logging.info("Finished extracting aligned sequences")
    
        # runs graftm for each of the HMMs doing the actual alignments, for each
        # of the input sequences
        logging.info("Running taxonomic assignment with graftm..")
        graftm_align_directory_base = os.path.join(working_directory, 'graftm_aligns')
        os.mkdir(graftm_align_directory_base)
        commands = []
        for sample_name in sample_names:
            if sample_name in sample_to_gpkg_to_input_sequences:
                for hmm_and_position in hmms:
                    key = hmm_and_position.gpkg_basename()
                    if key in sample_to_gpkg_to_input_sequences[sample_name]:
                        tmp_graft = sample_to_gpkg_to_input_sequences[sample_name][key]
                        cmd = "graftM graft --threads %i --verbosity 2 "\
                             "--forward %s "\
                             "--graftm_package %s --output_directory %s/%s_vs_%s "\
                             "--input_sequence_type nucleotide "\
                             "--assignment_method %s" % (\
                                    1, #use 1 thread since most likely better to parallelise processes with extern, not threads here
                                    tmp_graft.name,
                                    hmm_and_position.gpkg_path,
                                    graftm_align_directory_base,
                                    sample_name,
                                    hmm_and_position.gpkg_basename(),
                                    graftm_assignment_method)
                        if bootstrap_contigs:
                            bootstrap_hmm = bootstrap_hmms[hmm.hmm_filename]
                            if os.path.isfile(bootstrap_hmm):
                                cmd += " --search_hmm_files %s %s" % (
                                            bootstrap_hmm,
                                            hmm.hmm_path())
                        commands.append(cmd)
                    else:
                        logging.debug("No sequences found aligning from gpkg %s to sample %s, skipping" % (hmm_and_position.gpkg_basename(), sample_name)) 
            else:
                logging.debug("No sequences found aligning to sample %s at all, skipping" % sample_name)
        extern.run_many(commands, num_threads=num_threads)
        logging.info("Finished running taxonomic assignment with graftm")
    
        # get the sequences out for each of them
        with open(output_otu_table,'w') as output:
            to_print = split('gene sample sequence num_hits coverage taxonomy')
            if output_extras:
                to_print.append('read_names')
                to_print.append('nucleotides_aligned')
                if known_otu_tables:
                    to_print.append("taxonomy_by_known?")
            output.write("\t".join(to_print)+"\n")
            
            if known_otu_tables:
                logging.info("Parsing known taxonomy OTU tables")
                known_taxes = KnownOtuTable()
                known_taxes.parse_otu_tables(known_otu_tables)
                
            for sample_name in sample_names:
                if sample_name in sample_to_gpkg_to_input_sequences:
                    for hmm_and_position in hmms:
                        key = hmm_and_position.gpkg_basename()
                        if key in sample_to_gpkg_to_input_sequences[sample_name]:
                            tmp_graft = sample_to_gpkg_to_input_sequences[sample_name][key]
                            
                            tmpbase = os.path.basename(tmp_graft.name[:-6])#remove .fasta
                            base_dir = os.path.join(graftm_align_directory_base,
                                '%s_vs_%s' % (sample_name, hmm_and_position.gpkg_basename()),
                                tmpbase)
                            
                            proteins_file = os.path.join(base_dir, "%s_orf.fa" % tmpbase)
                            nucleotide_file = os.path.join(base_dir, "%s_hits.fa" % tmpbase)
                            aligned_seqs = self._get_windowed_sequences(
                                proteins_file,
                                nucleotide_file, hmm_and_position.hmm_path(),
                                hmm_and_position.best_position)
                            
                            if len(aligned_seqs) == 0:
                                logging.debug("Found no alignments for %s, skipping to next sample/hmm" % hmm_and_position.hmm_basename())
                                continue
                            logging.debug("Found %i sequences for hmm %s, sample '%s'" % (len(aligned_seqs),
                                                                                        hmm_and_position.hmm_basename(),
                                                                                        sample_name))
                            taxonomies = TaxonomyFile(os.path.join(base_dir, "%s_read_tax.tsv" % tmpbase))
            
                            # convert to OTU table, output
                            for info in self._seqs_to_counts_and_taxonomy(aligned_seqs,
                                                                    taxonomies):
                                if known_otu_tables:
                                    tax_assigned_through_known = False
                                to_print = [hmm_and_position.gpkg_basename(),
                                                sample_name,
                                                info.seq,
                                                str(info.count),
                                                "%.2f" % info.coverage]
                                if known_otu_tables and info.seq in known_taxes:
                                    tax_assigned_through_known = True
                                    to_print.append(known_taxes[info.seq].taxonomy)
                                else:
                                    to_print.append(info.taxonomy)
                                if output_extras:
                                    to_print.append(' '.join(info.names))
                                    to_print.append(' '.join([str(l) for l in info.aligned_lengths]))
                                    if known_otu_tables:
                                        to_print.append(tax_assigned_through_known)
                                output.write("\t".join(to_print) + "\n")
                    
        return_cleanly()
        
    def _get_windowed_sequences(self, protein_sequences_file, nucleotide_sequence_file, hmm_path, position):
        if os.stat(nucleotide_sequence_file).st_size == 0: return []
        nucleotide_sequences = SeqReader().read_nucleotide_sequences(nucleotide_sequence_file)
        protein_alignment = self._align_proteins_to_hmm(protein_sequences_file,
                                                      hmm_path)
        return MetagenomeOtuFinder().find_windowed_sequences(protein_alignment,
                                                        nucleotide_sequences,
                                                        20,
                                                        position)
        
    def _placement_input_fasta_name(self, hmm_and_position, sample_name, graftm_separate_directory_base):
        return '%s/%s_vs_%s/%s_hits/%s_hits_hits.fa' % (graftm_separate_directory_base,
                                                  sample_name,
                                                  os.path.basename(hmm_and_position.gpkg_path),
                                                  sample_name,
                                                  sample_name)
        
    def _extract_relevant_reads(self, graftm_separate_directory_base, sample_names, hmms):
        '''Given 'separates' directory, extract reads that will be used as
        part of the singlem choppage process as tempfiles in a hash'''
        sample_to_gpkg_to_input_sequences = {}
        for sample_name in sample_names:
            sample_to_gpkg_to_input_sequences[sample_name] = {}
            for hmm_and_position in hmms:
                base_dir = os.path.join(\
                    graftm_separate_directory_base,
                    "%s_vs_%s" % (sample_name,
                                  os.path.basename(hmm_and_position.gpkg_path)),
                    '%s_hits' % sample_name)
                protein_sequences_file = os.path.join(
                    base_dir, "%s_hits_orf.fa" % sample_name)
                nucleotide_sequence_fasta_file = os.path.join(
                    base_dir, "%s_hits_hits.fa" % sample_name)

                # extract the names of the relevant reads
                aligned_seqs = self._get_windowed_sequences(\
                        protein_sequences_file,
                        nucleotide_sequence_fasta_file,
                        hmm_and_position.hmm_path(),
                        hmm_and_position.best_position)
                if len(aligned_seqs) > 0:
                    tmp = tempfile.NamedTemporaryFile(prefix='singlem.%s.' % sample_name,suffix='.fasta')
                    cmd = "fxtract -X -H -f /dev/stdin %s > %s" % (nucleotide_sequence_fasta_file, tmp.name)
                    process = subprocess.Popen(['bash','-c',cmd],
                                               stdin=subprocess.PIPE)
                    process.communicate("\n".join([s.name for s in aligned_seqs]))
                    sample_to_gpkg_to_input_sequences[sample_name][os.path.basename(hmm_and_position.gpkg_path)] = tmp
        
        return sample_to_gpkg_to_input_sequences
        
    def _align_proteins_to_hmm(self, proteins_file, hmm_file):
        '''hmmalign proteins to hmm, and return an alignment object'''
        
        with tempfile.NamedTemporaryFile(prefix="singlem", suffix=".fasta") as f:
            cmd = "hmmalign %s %s |seqmagick convert --input-format stockholm - %s" % (hmm_file,
                                                              proteins_file,
                                                              f.name)
            logging.debug("Running cmd: %s" % cmd)
            extern.run(cmd)
            return SeqReader().protein_alignment_from_alignment_file(f.name)
    
    def _seqs_to_counts_and_taxonomy(self, sequences, taxonomies):
        '''given an array of Sequence objects, and hash of taxonomy file,
        yield over 'Info' objects that contain e.g. the counts of the aggregated
        sequences and corresponding median taxonomies.
        '''
        class CollectedInfo:
            def __init__(self):
                self.count = 0
                self.taxonomies = []
                self.names = []
                self.coverage = 0.0
                self.aligned_lengths = []
        
        seq_to_collected_info = {}
        for s in sequences:
            tax = taxonomies[s.name]
            try:
                collected_info = seq_to_collected_info[s.aligned_sequence]
            except KeyError:
                collected_info = CollectedInfo()
                seq_to_collected_info[s.aligned_sequence] = collected_info
                
            collected_info.count += 1
            if taxonomies: collected_info.taxonomies.append(tax)
            collected_info.names.append(s.name)
            collected_info.coverage += s.coverage_increment()
            collected_info.aligned_lengths.append(s.aligned_length)
        
        class Info:
            def __init__(self, seq, count, taxonomy, names, coverage, aligned_lengths):
                self.seq = seq
                self.count = count
                self.taxonomy = taxonomy
                self.names = names
                self.coverage = coverage
                self.aligned_lengths = aligned_lengths
            
        for seq, collected_info in seq_to_collected_info.iteritems():
            median_tax = self._median_taxonomy(collected_info.taxonomies)
            yield Info(seq,
                       collected_info.count,
                       median_tax,
                       collected_info.names,
                       collected_info.coverage,
                       collected_info.aligned_lengths)
    
    
    
    def _median_taxonomy(self, taxonomies):
        levels_to_counts = []
        for tax_string in taxonomies:
            for i, tax in enumerate(tax_string.split('; ')):
                if i >= len(levels_to_counts):
                    levels_to_counts.append({})
                try:
                    levels_to_counts[i][tax] += 1
                except KeyError:
                    levels_to_counts[i][tax] = 1
                    
                    
        median_tax = []
        for level_counts in levels_to_counts:
            max_count = 0
            max_tax = None
            for tax, count in level_counts.iteritems():
                if count > max_count:
                    max_count = count
                    max_tax = tax
            if float(max_count) / len(taxonomies) >= 0.5:
                median_tax.append(max_tax)
            else:
                break
        return '; '.join(median_tax)
                    