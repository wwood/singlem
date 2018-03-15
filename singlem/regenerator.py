import logging

from singlem.graftm_result import GraftMResult
from singlem.singlem_package import SingleMPackage

class Rgenerator:
    def regenerate(self, **kwargs):
        input_singlem_package = kwargs.pop('input_singlem_package')
        output_singlem_package = kwargs.pop('output_singlem_package')
        working_directory = kwargs.pop('working_directory')
        euk_sequences = kwargs.pop('euk_sequences')
        euk_taxonomy = kwargs.pop('euk_taxonomy')
        intermediate_archaea_graftm_package = kwargs.pop('intermediate_archaea_graftm_package')
        intermediate_bacteria_graftm_package = kwargs.pop('intermediate_bacteria_graftm_package')
        intermediate_sequences_for_inclusion = kwargs.pop('intermediate_sequences_for_inclusion')
        input_taxonomy = kwargs.pop('input_taxonomy')

        original_pkg = SingleMPackage.acquire(input_singlem_package)
        original_hmm_path = original_pkg.hmm_path()
        basename = original_pkg.graftm_package_basename()

        # Run GraftM on the euk sequences with the bacterial set
        euk_graftm_output = os.path.join(working_directory,
                                         "%s-euk_graftm" % basename)
        cmd = "graftM graft --graftm_package '%s' --sequences '%s' --output %s --force" % (
            original_pkg.graftm_package_path(),
            euk_sequences,
            euk_graftm_output)
        extern.run(cmd)

        # Extract hit sequences from that set
        euk_result = GraftMResult(euk_graftm_output)
        hit_paths = euk_result.unaligned_sequence_paths()
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
        with open(final_sequences_path, 'w') as final_seqs_fp:
            with open(euk_hits_path) as euk_seqs_fp:
                for name, seq, _ in SeqReader.readfq(euk_seqs_fp):
                    num_euk_hits += 1
                    #TODO: Dereplicate at some level
                    final_seqs_fp.write(">%s\n%s\n" % (name, seq))
        logging.info("Found %i eukaryotic sequences to include in the package" % \
                     num_euk_hits)
        archeal_seqs = archaeal_intermediate_pkg.unaligned_sequence_database_path()
        bacterial_seqs = bacterial_intermediate_pkg.unaligned_sequence_database_path()
        extern.run("cat %s >> %s" % (archeal_seqs, final_sequences_path))
        extern.run("cat %s >> %s" % (bacterial_seqs, final_sequences_path))

        # Concatenate euk and input taxonomy
        final_taxonomy_file = os.path.join(working_directory,
                                            "%s_final_taxonomy.csv" % basename)
        extern.run("cat %s %s > %s" % (
            euk_taxonomy, input_taxonomy, final_taxonomy_file))

        # Run graftm create to get the final package
        final_gpkg = os.path.join(working_directory,
                                  "%s_final.gpkg" % basename)
        extern.run("graftM create --sequences %s --taxonomy %s --search_hmm_files %s %s --hmm %s --output %s" % (
            final_sequences_path,
            final_taxonomy_file,
            ' '.join(archaeal_intermediate_pkg.search_hmm_paths),
            ' '.join(bacterial_intermediate_pkg.search_hmm_paths),
            original_hmm_path,
            final_gpkg))

        # Run singlem create to put the final package together
        SingleMPackage().compile(
            output_singlem_package,
            final_gpkg,
            original_pkg.singlem_position)


