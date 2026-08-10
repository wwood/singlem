import logging
import os
import csv
import glob
import re
import shutil

import extern

_GENOME_ACCESSION_REGEX = re.compile(r'(GC[AF]_\d+\.\d+)')

# weebill's two-stage databases, the only kind 'weebill profile --two-stage' reads
# and so the only kind a metapackage may bundle.
TWO_STAGE_DB_EXTENSION = '.syl2db'
# weebill writes compressed sample sketches.
SKETCH_EXTENSION = '.sylspc'

# weebill reports Eff_cov normally and renames it True_cov under -u
# (--estimate-unknown), which corrects for the unknown sequence fraction and so
# estimates the genome's actual coverage rather than one deflated by it. The name is
# the only signal of which was run, so it is preserved through annotation; condense
# reads either, since a profile supplied by hand to condense --weebill-profile need
# not have been made with -u.
UNKNOWN_CORRECTED_COVERAGE_COLUMN = 'True_cov'
COVERAGE_COLUMNS = [UNKNOWN_CORRECTED_COVERAGE_COLUMN, 'Eff_cov']


class WeebillProfiler:
    '''Runs weebill against the two-stage database(s) bundled in a metapackage and
    annotates its genome-level output with GTDB taxonomy, producing a TSV that
    condense's --weebill-profile path consumes (columns Sample_file, taxonomy, and
    the coverage column carried through under its original name).

    weebill is always run with -u (--estimate-unknown), under which it corrects each
    genome's coverage for the sample's unknown sequence fraction and reports it as
    True_cov rather than Eff_cov. That correction is what puts the coverages in
    SingleM's units, so condense can take alpha as 1 rather than calibrating it.

    Requires the weebill binary on the PATH.'''

    def run_from_reads(self, forward_reads, reverse_reads, metapackage, threads,
                       output_annotated_tsv, working_directory, sketch_output=None):
        '''Profile reads against the metapackage's database(s) and annotate. Returns
        output_annotated_tsv.

        Sketching is only done as a separate step when the sketches are wanted on
        disk (sketch_output), since otherwise profiling straight from the reads
        saves writing them and reading them straight back.'''
        databases, c = self._databases(metapackage)
        if sketch_output is None:
            raw_tsv = self._profile_reads(
                forward_reads, reverse_reads, databases, c, threads,
                os.path.join(working_directory, 'weebill_profile.tsv'))
        else:
            sketch_directory = os.path.join(working_directory, 'sketch')
            os.makedirs(sketch_directory, exist_ok=True)
            sketches = self._sketch_reads(forward_reads, reverse_reads, c, threads, sketch_directory)
            self._save_sketches(sketches, sketch_output)
            raw_tsv = self._profile_sketches(
                sketches, databases, threads,
                os.path.join(working_directory, 'weebill_profile.tsv'))
        return self._check_and_annotate(raw_tsv, metapackage, output_annotated_tsv)

    def run_from_sketch(self, sketch_path, metapackage, threads, output_annotated_tsv,
                        working_directory):
        '''Profile previously-saved sketch(es) against the metapackage's database(s)
        and annotate, so no access to the raw reads is needed.'''
        databases, _ = self._databases(metapackage)
        raw_tsv = self._profile_sketches(
            self._saved_sketches(sketch_path), databases, threads,
            os.path.join(working_directory, 'weebill_profile.tsv'))
        return self._check_and_annotate(raw_tsv, metapackage, output_annotated_tsv)

    def annotate(self, raw_profile_tsv, metapackage, output_tsv):
        '''Annotate a raw weebill profile with GTDB taxonomy from the metapackage,
        writing Sample_file / taxonomy / True_cov to output_tsv.'''
        rows = []
        needed = set()
        with open(raw_profile_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields = reader.fieldnames
            coverage_column = next((c for c in COVERAGE_COLUMNS if fields and c in fields), None)
            if fields is None or 'Genome_file' not in fields or coverage_column is None:
                raise Exception("Unexpected weebill profile format: {}".format(fields))
            for row in reader:
                accession = self._extract_accession(row['Genome_file'])
                rows.append((row.get('Sample_file', ''), accession, row[coverage_column], row['Genome_file']))
                if accession is not None:
                    needed.add(accession)

        accession_to_taxonomy = metapackage.genome_accession_to_taxonomy(needed)

        num_written = 0
        num_skipped = 0
        with open(output_tsv, 'w') as out:
            out.write('Sample_file\ttaxonomy\t{}\n'.format(coverage_column))
            for sample, accession, coverage, genome_file in rows:
                taxonomy = accession_to_taxonomy.get(accession)
                if taxonomy is None:
                    logging.debug("No metapackage taxonomy for weebill genome {} (accession {}), skipping".format(
                        genome_file, accession))
                    num_skipped += 1
                    continue
                out.write('\t'.join([sample, taxonomy, coverage]) + '\n')
                num_written += 1
        logging.info("Annotated {} weebill genome coverages with GTDB taxonomy ({} unmatched)".format(
            num_written, num_skipped))
        return output_tsv

    def _databases(self, metapackage):
        '''The bundled databases and the single -c they were all built at.'''
        databases = metapackage.weebill_databases()
        if len(databases) == 0:
            raise Exception("Metapackage bundles no weebill databases")
        return [db for db, _ in databases], databases[0][1]

    def _sketch_reads(self, forward_reads, reverse_reads, c, threads, output_directory):
        '''Sketch reads into output_directory; return the sorted list of sketch paths
        produced (one per sample).'''
        if not forward_reads:
            raise Exception("No reads provided to weebill sketch")
        logging.info("Sketching reads with weebill ..")
        extern.run("weebill sketch -c {} -t {} --compressed-database {}{}".format(
            c, threads, output_directory, self._read_arguments(forward_reads, reverse_reads)))
        sketches = sorted(glob.glob(os.path.join(output_directory, '*' + SKETCH_EXTENSION)))
        if len(sketches) == 0:
            raise Exception("weebill sketch produced no {} file".format(SKETCH_EXTENSION))
        return sketches

    def _profile_reads(self, forward_reads, reverse_reads, databases, c, threads, output_tsv):
        '''Profile reads directly, letting weebill sketch them at -c as it goes.'''
        if not forward_reads:
            raise Exception("No reads provided to weebill profile")
        logging.info("Profiling reads against {} weebill database(s) ..".format(len(databases)))
        extern.run("weebill profile --two-stage {} -u -c {} -t {} -o {}{}".format(
            ' '.join(databases), c, threads, output_tsv,
            self._read_arguments(forward_reads, reverse_reads)))
        return output_tsv

    def _profile_sketches(self, sketches, databases, threads, output_tsv):
        '''Profile sample sketches. -c is baked into them, so it is not repeated.
        Passing the databases together lets weebill reassign shared k-mers across
        them.'''
        logging.info("Profiling {} sample sketch(es) against {} weebill database(s) ..".format(
            len(sketches), len(databases)))
        extern.run("weebill profile --two-stage -u -t {} -o {} {} {}".format(
            threads, output_tsv, ' '.join(databases), ' '.join(sketches)))
        return output_tsv

    def _read_arguments(self, forward_reads, reverse_reads):
        '''The -1/-2 (or -r) fragment naming the read files, shared by sketching and
        by profiling straight from reads.'''
        if reverse_reads:
            return " -1 {} -2 {}".format(' '.join(forward_reads), ' '.join(reverse_reads))
        return " -r {}".format(' '.join(forward_reads))

    def _save_sketches(self, sketches, sketch_output):
        os.makedirs(sketch_output, exist_ok=True)
        for s in sketches:
            shutil.copy(s, os.path.join(sketch_output, os.path.basename(s)))
        logging.info("Saved {} weebill sketch(es) to {}/".format(len(sketches), sketch_output))

    def _saved_sketches(self, sketch_path):
        '''The sketch(es) at sketch_path, which is either a directory written by
        _save_sketches or a single sketch file.'''
        if os.path.isdir(sketch_path):
            sketches = sorted(glob.glob(os.path.join(sketch_path, '*' + SKETCH_EXTENSION)))
        else:
            sketches = [sketch_path]
        if len(sketches) == 0:
            raise Exception("No weebill sketch ({}) found at {}".format(SKETCH_EXTENSION, sketch_path))
        return sketches

    def _check_and_annotate(self, raw_profile_tsv, metapackage, output_annotated_tsv):
        '''Fail loudly if -u did not produce the unknown-corrected coverage it was
        asked for. Callers take alpha as 1 on the strength of that correction, so a
        silently uncorrected profile would be read on the wrong scale.'''
        with open(raw_profile_tsv) as f:
            fields = csv.DictReader(f, delimiter='\t').fieldnames or []
        if UNKNOWN_CORRECTED_COVERAGE_COLUMN not in fields:
            raise Exception(
                "weebill was run with -u (--estimate-unknown) but its profile reports no {} "
                "column (found {}). Its coverages are then not on SingleM's scale and cannot "
                "be used with alpha=1.".format(UNKNOWN_CORRECTED_COVERAGE_COLUMN, fields))
        return self.annotate(raw_profile_tsv, metapackage, output_annotated_tsv)

    def _extract_accession(self, genome_file):
        match = _GENOME_ACCESSION_REGEX.search(genome_file)
        return match.group(1) if match else None
