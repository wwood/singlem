import logging
import os

import extern
import dendropy
from graftm.graftm_package import GraftMPackage

from graftm_result import GraftMResult
from singlem_package import SingleMPackageVersion2, SingleMPackage
from sequence_classes import SeqReader
from dereplicator import Dereplicator
from sequence_extractor import SequenceExtractor

class Regenerator:
    def regenerate(self, **kwargs):
        input_singlem_package = kwargs.pop('input_singlem_package')
        output_singlem_package = kwargs.pop('output_singlem_package')
        working_directory = kwargs.pop('working_directory')
        euk_sequences = kwargs.pop('euk_sequences')
        euk_taxonomy = kwargs.pop('euk_taxonomy')
        intermediate_archaea_graftm_package = kwargs.pop('intermediate_archaea_graftm_package')
        intermediate_bacteria_graftm_package = kwargs.pop('intermediate_bacteria_graftm_package')
        input_taxonomy = kwargs.pop('input_taxonomy')
        type_strains_list_file = kwargs.pop('type_strains_list_file')

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        original_pkg = SingleMPackage.acquire(input_singlem_package)
        original_hmm_path = original_pkg.hmm_path()
        basename = original_pkg.graftm_package_basename()

        # Run GraftM on the euk sequences with the bacterial set
        euk_graftm_output = os.path.join(working_directory,
                                         "%s-euk_graftm" % basename)
        cmd = "graftM graft --graftm_package '%s' --search_and_align_only --forward '%s' --output %s --force" % (
            original_pkg.graftm_package_path(),
            euk_sequences,
            euk_graftm_output)
        extern.run(cmd)

        # Extract hit sequences from that set
        euk_result = GraftMResult(euk_graftm_output)
        hit_paths = euk_result.unaligned_sequence_paths(require_hits=True)
        if len(hit_paths) != 1: raise Exception(
                "Unexpected number of hits against euk in graftm")
        euk_hits_path = hit_paths.values()[0]

        # Concatenate euk, archaea and bacterial sequences
        archaeal_intermediate_pkg = GraftMPackage.acquire(
            intermediate_archaea_graftm_package)
        bacterial_intermediate_pkg = GraftMPackage.acquire(
            intermediate_bacteria_graftm_package)
        num_euk_hits = 0
        final_sequences_path = os.path.join(working_directory,
                                            "%s_final_sequences.faa" % basename)
        archeal_seqs = archaeal_intermediate_pkg.unaligned_sequence_database_path()
        bacterial_seqs = bacterial_intermediate_pkg.unaligned_sequence_database_path()
        with open(type_strains_list_file) as f:
            type_strain_identifiers = [s.strip() for s in f.readlines()]
        logging.info("Read in %i type strain IDs e.g. %s" % (
            len(type_strain_identifiers), type_strain_identifiers[0]))

        with open(final_sequences_path, 'w') as final_seqs_fp:
            with open(euk_hits_path) as euk_seqs_fp:
                for name, seq, _ in SeqReader().readfq(euk_seqs_fp):
                    if name.find('_split_') == -1:
                        num_euk_hits += 1
                        #TODO: Dereplicate at some level
                        final_seqs_fp.write(">%s\n%s\n" % (name, seq))
            logging.info("Found %i eukaryotic sequences to include in the package" % \
                         num_euk_hits)

            # Dereplicate hit sequences on the species level, choosing type strains
            # where applicable.
            dereplicator = Dereplicator()
            for gpkg in [archaeal_intermediate_pkg, bacterial_intermediate_pkg]:
                tax = gpkg.taxonomy_hash()
                species_dereplicated_ids = dereplicator.dereplicate(
                    list(tax.keys()),
                    8, # root, kingdom, phylum, c o f g s
                    tax,
                    type_strain_identifiers)
                logging.debug("Dereplicator returned %i entries" % len(species_dereplicated_ids))
                num_total = 0
                num_written = 0
                with open(gpkg.unaligned_sequence_database_path()) as seqs:
                    for name, seq, _ in SeqReader().readfq(seqs):
                        num_total += 1
                        if name in species_dereplicated_ids:
                            final_seqs_fp.write(">%s\n%s\n" % (name, seq))
                            num_written += 1
                logging.info(
                    "Of %i sequences in gpkg %s, %i species-dereplicated were included in the final package." %(
                        num_total, gpkg, num_written))

        # Concatenate euk and input taxonomy
        final_taxonomy_file = os.path.join(working_directory,
                                            "%s_final_taxonomy.csv" % basename)
        extern.run("cat %s %s > %s" % (
            euk_taxonomy, input_taxonomy, final_taxonomy_file))

        # Run graftm create to get the final package
        final_gpkg = os.path.join(working_directory,
                                  "%s_final.gpkg" % basename)
        cmd = "graftM create --force --sequences %s --taxonomy %s --search_hmm_files %s %s --hmm %s --output %s" % (
            final_sequences_path,
            final_taxonomy_file,
            ' '.join(archaeal_intermediate_pkg.search_hmm_paths()),
            ' '.join(bacterial_intermediate_pkg.search_hmm_paths()),
            original_hmm_path,
            final_gpkg)
        extern.run(cmd)

        ##############################################################################
        # Remove sequences from the diamond DB that are not in the tree i.e.
        # those that are exact duplicates, so that the diamond_example hits are
        # always in the tree.
        # Read the list of IDs in the tree with dendropy
        final_gpkg_object = GraftMPackage.acquire(final_gpkg)
        unaligned_seqs = final_gpkg_object.unaligned_sequence_database_path()
        tree = dendropy.Tree.get(path=final_gpkg_object.reference_package_tree_path(),
                                 schema='newick')
        leaf_names = [l.taxon.label.replace(' ','_') for l in tree.leaf_node_iter()]
        logging.debug("Read in final tree with %i leaves" % len(leaf_names))

        # Extract out of the sequences file in the graftm package
        final_seqs = SequenceExtractor().extract_and_read(
            leaf_names, unaligned_seqs)
        if len(final_seqs) != len(leaf_names):
            raise Exception("Do not appear to have extracted the expected number of sequences from the unaligned fastat file")

        # Write the reads into sequences file in place
        with open(unaligned_seqs, 'w') as f:
            for s in final_seqs:
                f.write(">%s\n" % s.name)
                f.write(s.seq)
                f.write("\n")

        # Regenerate the diamond DB
        final_gpkg_object.create_diamond_db()

        ##############################################################################
        # Run singlem create to put the final package together
        SingleMPackageVersion2.compile(
            output_singlem_package,
            final_gpkg,
            original_pkg.singlem_position(),
            original_pkg.window_size())
        logging.info("SingleM package generated.")


