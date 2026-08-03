import logging
import os

import extern

from .sylph import SylphProfiler

# weebill is a fork of sylph. Its two-stage databases carry a distinct extension,
# and are the only kind 'weebill profile --two-stage' accepts, so the extension is
# what tells a bundled database's profiler apart from a plain sylph .syldb.
TWO_STAGE_DB_EXTENSION = '.syl2db'
# weebill writes compressed sample sketches (.sylspc) where sylph writes .sylsp.
WEEBILL_SKETCH_EXTENSION = '.sylspc'


def is_two_stage_database(path):
    '''True if path names a weebill two-stage database rather than a sylph one.'''
    return path.endswith(TWO_STAGE_DB_EXTENSION)


def read_profiler_for_metapackage(metapackage):
    '''Return the profiler for the databases the metapackage bundles -- weebill for
    two-stage (.syl2db) databases, sylph for plain .syldb -- or None when it bundles
    none. Raises if the two kinds are mixed, since neither binary reads both.'''
    databases = metapackage.sylph_databases()
    if len(databases) == 0:
        return None
    two_stage = [is_two_stage_database(db) for db, _ in databases]
    if all(two_stage):
        return WeebillProfiler()
    if any(two_stage):
        raise Exception(
            "Metapackage bundles a mix of weebill two-stage ({}) and sylph databases, which "
            "cannot be profiled together: {}".format(
                TWO_STAGE_DB_EXTENSION, [db for db, _ in databases]))
    return SylphProfiler()


class WeebillProfiler(SylphProfiler):
    '''Runs weebill against the two-stage (.syl2db) database(s) bundled in a
    metapackage and annotates its genome-level output with GTDB taxonomy, producing
    the same TSV that condense's --sylph-profile path consumes.

    weebill is a fork of sylph, so its profile output is sylph's and annotation is
    inherited unchanged. What differs is how it is run: 'profile --two-stage' reads
    the seekable database rather than loading it, and it will sketch reads itself,
    so the common case (reads in, profile out) is one invocation rather than a
    sketch followed by a profile.

    Requires the weebill binary on the PATH.'''

    BINARY = 'weebill'
    SKETCH_EXTENSION = WEEBILL_SKETCH_EXTENSION

    def sketch_reads(self, forward_reads, reverse_reads, c, threads, output_directory):
        '''Sketch reads into output_directory, in the compressed form the two-stage
        profiler reads; return the sorted list of sketch paths produced.'''
        if not forward_reads:
            raise Exception("No reads provided to weebill sketch")
        cmd = "weebill sketch -c {} -t {} --compressed-database {}".format(
            c, threads, output_directory)
        cmd += self._read_arguments(forward_reads, reverse_reads)
        logging.info("Sketching reads with weebill ..")
        extern.run(cmd)
        return self._sketches_in_directory(output_directory)

    def profile(self, sketch_paths, sylph_dbs, threads, output_tsv, estimate_unknown=False):
        '''Profile sample sketches against one or more two-stage databases sharing
        the -c baked into the sketches.'''
        logging.info("Profiling {} sample sketch(es) against {} weebill database(s) ..".format(
            len(sketch_paths), len(sylph_dbs)))
        extern.run("weebill profile --two-stage{} -t {} -o {} {} {}".format(
            ' -u' if estimate_unknown else '', threads, output_tsv,
            ' '.join(sylph_dbs), ' '.join(sketch_paths)))
        return output_tsv

    def profile_reads(self, forward_reads, reverse_reads, sylph_dbs, c, threads, output_tsv,
                      estimate_unknown=False):
        '''Profile reads directly, letting weebill sketch them at -c as it goes.'''
        if not forward_reads:
            raise Exception("No reads provided to weebill profile")
        logging.info("Profiling reads against {} weebill database(s) ..".format(len(sylph_dbs)))
        cmd = "weebill profile --two-stage{} -c {} -t {} -o {}".format(
            ' -u' if estimate_unknown else '', c, threads, output_tsv)
        cmd += self._read_arguments(forward_reads, reverse_reads)
        cmd += " {}".format(' '.join(sylph_dbs))
        extern.run(cmd)
        return output_tsv

    def run_from_reads(self, forward_reads, reverse_reads, metapackage, threads,
                       output_annotated_tsv, working_directory, sketch_output=None,
                       estimate_unknown=True):
        '''Profile reads against the metapackage's two-stage database(s) and annotate.

        Sketching is only done as a separate step when the sketches are wanted on
        disk (sketch_output), since otherwise profiling straight from the reads
        saves writing and re-reading them.'''
        if sketch_output is not None:
            return super().run_from_reads(
                forward_reads, reverse_reads, metapackage, threads, output_annotated_tsv,
                working_directory, sketch_output=sketch_output, estimate_unknown=estimate_unknown)

        databases_by_c = self._group_databases_by_c(self._databases(metapackage))
        raw_tsvs = []
        for i, c in enumerate(sorted(databases_by_c)):
            raw_tsv = os.path.join(working_directory, 'weebill_profile_{}.tsv'.format(i))
            self.profile_reads(forward_reads, reverse_reads, databases_by_c[c], c, threads,
                               raw_tsv, estimate_unknown=estimate_unknown)
            raw_tsvs.append(raw_tsv)
        merged = self._merge_raw_profiles(
            raw_tsvs, os.path.join(working_directory, 'weebill_profile_merged.tsv'))
        if estimate_unknown:
            self._check_unknown_corrected(merged)
        return self.annotate(merged, metapackage, output_annotated_tsv)

    def run_from_sketch(self, sketch_path, metapackage, threads, output_annotated_tsv,
                        working_directory, estimate_unknown=True):
        return super().run_from_sketch(
            sketch_path, metapackage, threads, output_annotated_tsv, working_directory,
            estimate_unknown=estimate_unknown)

    def _databases(self, metapackage):
        databases = metapackage.sylph_databases()
        if len(databases) == 0:
            raise Exception("Metapackage bundles no weebill databases")
        return databases
