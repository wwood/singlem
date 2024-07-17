import logging
import os
import tempfile
import extern
import pickle
import shutil

from graftm.graftm_package import GraftMPackage

from .graftm_result import GraftMResult
from .singlem_package import SingleMPackageVersion4, SingleMPackage
from .sequence_classes import SeqReader, Sequence
from .metagenome_otu_finder import MetagenomeOtuFinder
from .pipe_sequence_extractor import _align_proteins_to_hmm
from . import CREATE_MIN_ALIGNED_PERCENT

class Regenerator:

    def regenerate(self, **kwargs):
        input_singlem_package = kwargs.pop('input_singlem_package')
        output_singlem_package = kwargs.pop('output_singlem_package')
        input_sequences = kwargs.pop('input_sequences')
        input_taxonomy = kwargs.pop('input_taxonomy')
        euk_sequences = kwargs.pop('euk_sequences', None) # This is called euk_sequences for historical reasons, but is actually just any off-target seqs
        euk_taxonomy = kwargs.pop('euk_taxonomy', None) # Likewise
        window_position = kwargs.pop('window_position', None)
        sequence_prefix = kwargs.pop('sequence_prefix', None)
        min_aligned_percent = kwargs.pop('min_aligned_percent', CREATE_MIN_ALIGNED_PERCENT)
        no_further_euks = kwargs.pop('no_further_euks', None)

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        original_pkg = SingleMPackage.acquire(input_singlem_package)
        original_hmm_path = original_pkg.hmm_path()
        original_hmm_search_paths = original_pkg.graftm_package().search_hmm_paths()
        basename = original_pkg.graftm_package_basename()

        if window_position:
            output_window_position = window_position
        else:
            output_window_position = original_pkg.singlem_position()

        # Ensure protein package
        if not original_pkg.is_protein_package():
            raise Exception("Only works with protein packages.")

        # Ensure v3 or above
        if not original_pkg.version >= 3:
            raise Exception("Only works with v3 and above packages for the moment.")

        tmp = tempfile.TemporaryDirectory()
        working_directory = tmp.name
        logging.debug("Using working directory %s" % working_directory)

        # Files to concatenate off_target and input seqs/taxonomy
        final_taxonomy_path = os.path.join(working_directory, "%s_final_taxonomy.csv" % basename)
        final_sequences_path = os.path.join(working_directory, "%s_final_sequences.faa" % basename)

        if no_further_euks:
            logging.info("Skipping off-target search. Copying input files directly")
            shutil.copyfile(input_taxonomy, final_taxonomy_path)
            shutil.copyfile(input_sequences, final_sequences_path)

        else:
            # Run GraftM on the euk sequences with the original graftm package
            logging.info("Finding off-target sequences via GraftM ..")
            euk_graftm_output = os.path.join(working_directory, "%s-euk_graftm" % basename)
            cmd = "graftM graft --graftm_package '%s' --search_and_align_only --forward '%s' --output %s --force" % (
                original_pkg.graftm_package_path(), euk_sequences, euk_graftm_output)
            extern.run(cmd)

            # Extract hit sequences from that set
            euk_result = GraftMResult(euk_graftm_output, False)
            hit_paths = euk_result.unaligned_sequence_paths(require_hits=True)

            with open(final_taxonomy_path, 'w') as final_tax_fp:
                # Copy across input taxonomy
                with open(input_taxonomy) as input_tax_f:
                    final_tax_fp.write(input_tax_f.read())
                    
                    if len(hit_paths) == 0:
                        euk_taxonomy = ""
                        logging.info("Found no off-target sequences that match")
                    else:
                        euk_hits_path = next(iter(hit_paths.values()))  #i.e. first

                        # Concatenate input (bacteria and archaea) and euk sequences
                        num_euk_hits = 0

                        # Read in taxonomy hash
                        euk_name_to_taxonomy = {}
                        with open(euk_taxonomy) as f:
                            for line in f:
                                splits = line.split('\t')
                                if len(splits) != 2:
                                    raise Exception("Unexpected taxonomy line found: {}".format(line))
                                euk_name_to_taxonomy[splits[0].strip()] = splits[1].strip()

                            with open(final_sequences_path, 'w') as final_seqs_fp:
                                with open(euk_hits_path) as euk_seqs_fp:
                                    for name, seq, _ in SeqReader().readfq(euk_seqs_fp):
                                        if name.find('_split_') == -1:
                                            num_euk_hits += 1
                                            final_seqs_fp.write(">%s\n%s\n" % (name, seq))
                                            final_tax_fp.write("{}\t{}\n".format(name, euk_name_to_taxonomy[name]))
                        logging.info("Found %i off-target sequences to include in the package" % num_euk_hits)

            extern.run("cat %s >> %s" % (input_sequences, final_sequences_path))

        # Add package prefix to sequences and taxonomy
        if sequence_prefix != "" and sequence_prefix is not None:
            extern.run("sed -i 's/>/>{}/g' {}".format(sequence_prefix, final_sequences_path))
            extern.run("sed -i 's/^/{}/g' {}".format(sequence_prefix, final_taxonomy_path))

        # Run graftm create to get the output package
        final_gpkg = os.path.join(working_directory, "%s_final.gpkg" % basename)
        cmd = "graftM create --no-tree --force --min_aligned_percent %s --sequences %s --taxonomy %s --search_hmm_files %s --hmm %s --output %s" % (
            min_aligned_percent, final_sequences_path, final_taxonomy_path, ' '.join(original_hmm_search_paths),
            original_hmm_path, final_gpkg)

        extern.run(cmd)
        output_gpkg = GraftMPackage.acquire(final_gpkg)

        # Trim unaligned sequences according to alignment window +/- 30 aa
        logging.info("Trimming unaligned sequences according to alignment window")
        unaligned_basename = os.path.basename(output_gpkg.unaligned_sequence_database_path())
        protein_sequences = SeqReader().readfq(open(os.path.join(final_gpkg, unaligned_basename)))
        tmp_alignment = _align_proteins_to_hmm(protein_sequences, output_gpkg.alignment_hmm_path())
        logging.info("Found {} total sequences in alignment".format(len(tmp_alignment)))
        best_position = original_pkg.singlem_position()
        stretch_length = original_pkg.window_size() / 3
        trimmed_output = []

        ignored_columns = MetagenomeOtuFinder()._find_lower_case_columns(tmp_alignment)
        logging.debug("Ignoring columns %s", str(ignored_columns))
        # Find start of window in aligned sequence
        start_position = MetagenomeOtuFinder()._upper_case_position_to_alignment_position(
            best_position, ignored_columns)
        logging.debug("Using pre-defined best section of the alignment starting from %i" % (start_position + 1))
        # Find all positions in the window
        chosen_positions = MetagenomeOtuFinder()._best_position_to_chosen_positions(start_position, stretch_length,
                                                                                    ignored_columns)
        logging.debug("Found chosen positions %s", chosen_positions)

        for aligned_sequence in tmp_alignment:
            windowed_residues = aligned_sequence.seq[min(chosen_positions):1 + max(chosen_positions)].replace('-', '')
            residues_before_window = aligned_sequence.seq[0:min(chosen_positions)].replace('-', '')
            residues_after_window = aligned_sequence.seq[(max(chosen_positions) + 1):].replace('-', '')

            prepending_index = len(residues_before_window) - 30
            if prepending_index < 0:
                prepending_index = 0
            prepending_residues = residues_before_window[prepending_index:]
            final_index = 30
            if len(residues_after_window) < 30:
                final_index = len(residues_after_window)
            appending_residues = residues_after_window[0:final_index]

            final_sequence = prepending_residues + windowed_residues + appending_residues

            trimmed_output.append(Sequence(aligned_sequence.name, final_sequence.upper()))

        unaligned_graftm_file = os.path.join(final_gpkg, unaligned_basename)
        with open(unaligned_graftm_file, "w") as out:
            for entry in trimmed_output:
                out.write(entry.fasta())

        # Recreate diamond DB based on the trimmed sequences and add makeidx
        logging.info("Recreating diamond database")
        cmd = "diamond makedb --in %s --db %s" % (unaligned_graftm_file, output_gpkg.diamond_database_path())
        extern.run(cmd)
        logging.info("Adding makeidx to diamond database")
        cmd = "diamond makeidx -d %s" % output_gpkg.diamond_database_path()
        extern.run(cmd)

        # Create taxonomy hash
        logging.debug("Creating taxonomy hash pickle")
        taxonomy_hash_path = os.path.join(working_directory, "taxonomy_hash.pickle")
        tax_hash = {}
        with open(taxonomy_hash_path, 'wb') as file:
            with open(final_taxonomy_path) as tax:
                for line in tax:
                    split_line = line.strip().split('\t')
                    tax_hash[split_line[0]] = [taxa.strip() for taxa in split_line[1].split(';')]

            pickle.dump(tax_hash, file)

        ##############################################################################
        # Run singlem create to put the final package together
        SingleMPackageVersion4.compile(output_singlem_package, final_gpkg, output_window_position,
                                       original_pkg.window_size(), original_pkg.target_domains(),
                                       original_pkg.gene_description(), taxonomy_hash_path)
        logging.info("SingleM package generated.")
