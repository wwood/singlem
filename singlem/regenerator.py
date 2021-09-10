import logging
import os
import tempfile
import extern
import dendropy
from graftm.graftm_package import GraftMPackage

from .graftm_result import GraftMResult
from .singlem_package import SingleMPackageVersion3, SingleMPackageVersion2, SingleMPackage
from .sequence_classes import SeqReader
from .dereplicator import Dereplicator
from .sequence_extractor import SequenceExtractor


class Regenerator:
    def regenerate(self, **kwargs):
        input_singlem_package = kwargs.pop('input_singlem_package')
        output_singlem_package = kwargs.pop('output_singlem_package')
        input_sequences = kwargs.pop('input_sequences')
        input_taxonomy = kwargs.pop('input_taxonomy')
        euk_sequences = kwargs.pop('euk_sequences')
        euk_taxonomy = kwargs.pop('euk_taxonomy')
        min_aligned_percent = kwargs.pop('min_aligned_percent')

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        original_pkg = SingleMPackage.acquire(input_singlem_package)
        original_hmm_path = original_pkg.hmm_path()
        original_hmm_search_paths = original_pkg.graftm_package().search_hmm_paths()
        basename = original_pkg.graftm_package_basename()

        tmp = tempfile.TemporaryDirectory()
        working_directory = tmp.name
        logging.debug("Using working directory %s" % working_directory)


        # Run GraftM on the euk sequences with the bacterial set
        euk_graftm_output = os.path.join(working_directory,
                                         "%s-euk_graftm" % basename)
        cmd = "graftM graft --graftm_package '%s' --search_and_align_only --forward '%s' --output %s --force" % (
            original_pkg.graftm_package_path(),
            euk_sequences,
            euk_graftm_output)
        extern.run(cmd)

        # Extract hit sequences from that set
        euk_result = GraftMResult(euk_graftm_output, False)
        hit_paths = euk_result.unaligned_sequence_paths(require_hits=True)
        if len(hit_paths) != 1: raise Exception(
                "Unexpected number of hits against euk in graftm")
        euk_hits_path = next(iter(hit_paths.values())) #i.e. first

        # Concatenate input (bacteria and archaea) and euk sequences
        num_euk_hits = 0
        final_sequences_path = os.path.join(working_directory,
                                            "%s_final_sequences.faa" % basename)

        with open(final_sequences_path, 'w') as final_seqs_fp:
            with open(euk_hits_path) as euk_seqs_fp:
                for name, seq, _ in SeqReader().readfq(euk_seqs_fp):
                    if name.find('_split_') == -1:
                        num_euk_hits += 1
                        final_seqs_fp.write(">%s\n%s\n" % (name, seq))
            logging.info("Found %i eukaryotic sequences to include in the package" % \
                         num_euk_hits)
        
        extern.run("cat %s >> %s" % (input_sequences, final_sequences_path))

        # Concatenate euk and input taxonomy
        final_taxonomy_file = os.path.join(working_directory,
                                            "%s_final_taxonomy.csv" % basename)
        extern.run("cat %s %s > %s" % (
            euk_taxonomy, input_taxonomy, final_taxonomy_file))

        # Run graftm create to get the final package
        final_gpkg = os.path.join(working_directory,
                                  "%s_final.gpkg" % basename)
        cmd = "graftM create --no-tree --force --min_aligned_percent %s --sequences %s --taxonomy %s --search_hmm_files %s --hmm %s --output %s" % (
            min_aligned_percent,
            final_sequences_path,
            final_taxonomy_file,
            ' '.join(original_hmm_search_paths),
            original_hmm_path,
            final_gpkg)

        extern.run(cmd)

        ##############################################################################
        # Run singlem create to put the final package together
        if original_pkg.version == 2:
            SingleMPackageVersion2.compile(
                output_singlem_package,
                final_gpkg,
                original_pkg.singlem_position(),
                original_pkg.window_size())
        elif original_pkg.version == 3:
            SingleMPackageVersion3.compile(
                output_singlem_package,
                final_gpkg,
                original_pkg.singlem_position(),
                original_pkg.window_size(),
                original_pkg.target_domains(),
                original_pkg.gene_description())
        logging.info("SingleM package generated.")


