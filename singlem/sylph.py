import logging
import os
import csv
import glob
import re
import shutil
import tempfile

import extern

_GENOME_ACCESSION_REGEX = re.compile(r'(GC[AF]_\d+\.\d+)')


class SylphProfiler:
    '''Runs sylph (sketch + profile) and annotates its genome-level output with
    GTDB taxonomy from the metapackage, producing a TSV that condense's
    --sylph-profile path consumes (columns Sample_file, taxonomy, and the
    coverage column carried through under its original name).

    Requires the sylph binary on the PATH.'''

    # The binary and sample sketch extension, overridden by WeebillProfiler,
    # which is a sylph fork with its own name and compressed sketch format.
    BINARY = 'sylph'
    SKETCH_EXTENSION = '.sylsp'

    # sylph reports Eff_cov normally and renames it True_cov under -u
    # (--estimate-unknown), which corrects for the unknown sequence fraction and
    # so estimates the genome's actual coverage rather than one deflated by it.
    # The name is the only signal of which was run, so it is preserved through
    # annotation; condense reads either.
    COVERAGE_COLUMNS = ['True_cov', 'Eff_cov']

    def sketch_reads(self, forward_reads, reverse_reads, c, threads, output_directory):
        '''Sketch reads into output_directory; return the sorted list of sketch
        paths produced (one per sample). reverse_reads may be None for single-end.'''
        if not forward_reads:
            raise Exception("No reads provided to {} sketch".format(self.BINARY))
        cmd = "{} sketch -c {} -t {} -d {}".format(self.BINARY, c, threads, output_directory)
        cmd += self._read_arguments(forward_reads, reverse_reads)
        logging.info("Sketching reads with {} ..".format(self.BINARY))
        extern.run(cmd)
        return self._sketches_in_directory(output_directory)

    def _read_arguments(self, forward_reads, reverse_reads):
        '''The -1/-2 (or -r) fragment naming the read files, shared by sketch and
        by profiling straight from reads.'''
        if reverse_reads:
            return " -1 {} -2 {}".format(' '.join(forward_reads), ' '.join(reverse_reads))
        return " -r {}".format(' '.join(forward_reads))

    def _sketches_in_directory(self, directory):
        sketches = sorted(glob.glob(os.path.join(directory, '*' + self.SKETCH_EXTENSION)))
        if len(sketches) == 0:
            raise Exception("{} sketch produced no {} file".format(self.BINARY, self.SKETCH_EXTENSION))
        return sketches

    def profile(self, sketch_paths, sylph_dbs, threads, output_tsv, estimate_unknown=False):
        '''Profile sample sketches against one or more sylph databases (all of
        which must share the -c baked into the sketches), writing the raw sylph
        TSV to output_tsv. Passing the databases together lets sylph reassign
        shared k-mers across them.'''
        logging.info("Profiling {} sample sketch(es) against {} {} database(s) ..".format(
            len(sketch_paths), len(sylph_dbs), self.BINARY))
        extern.run("{} profile{} {} {} -t {} -o {}".format(
            self.BINARY, ' -u' if estimate_unknown else '',
            ' '.join(sylph_dbs), ' '.join(sketch_paths), threads, output_tsv))
        return output_tsv

    def annotate(self, raw_profile_tsv, metapackage, output_tsv):
        '''Annotate a raw sylph profile with GTDB taxonomy from the metapackage,
        writing Sample_file / taxonomy / Eff_cov to output_tsv.'''
        rows = []
        needed = set()
        with open(raw_profile_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields = reader.fieldnames
            coverage_column = next((c for c in self.COVERAGE_COLUMNS if fields and c in fields), None)
            if fields is None or 'Genome_file' not in fields or coverage_column is None:
                raise Exception("Unexpected sylph profile format: {}".format(fields))
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
            for sample, accession, eff_cov, genome_file in rows:
                taxonomy = accession_to_taxonomy.get(accession)
                if taxonomy is None:
                    logging.debug("No metapackage taxonomy for sylph genome {} (accession {}), skipping".format(
                        genome_file, accession))
                    num_skipped += 1
                    continue
                out.write('\t'.join([sample, taxonomy, eff_cov]) + '\n')
                num_written += 1
        logging.info("Annotated {} sylph genome coverages with GTDB taxonomy ({} unmatched)".format(
            num_written, num_skipped))
        return output_tsv

    def run_from_reads(self, forward_reads, reverse_reads, metapackage, threads,
                       output_annotated_tsv, working_directory, sketch_output=None,
                       estimate_unknown=False):
        '''Sketch reads (once per distinct -c across the metapackage's sylph
        databases), optionally save the sketches, profile each database against
        the sketch made at its -c, merge, and annotate. Returns output_annotated_tsv.'''
        databases = metapackage.sylph_databases()
        if len(databases) == 0:
            raise Exception("Metapackage bundles no {} databases".format(self.BINARY))

        # Sketch reads once per distinct -c value.
        databases_by_c = self._group_databases_by_c(databases)
        sketches_by_c = {}
        for c in sorted(databases_by_c):
            sketch_dir = os.path.join(working_directory, 'sketch_c{}'.format(c))
            os.makedirs(sketch_dir, exist_ok=True)
            sketches_by_c[c] = self.sketch_reads(forward_reads, reverse_reads, c, threads, sketch_dir)
        if sketch_output is not None:
            self._save_sketches_by_c(sketches_by_c, sketch_output)

        # Profile all databases that share a -c together, so sylph can reassign
        # shared k-mers across them.
        raw_tsvs = []
        for i, c in enumerate(sorted(databases_by_c)):
            raw_tsv = os.path.join(working_directory, 'sylph_profile_{}.tsv'.format(i))
            self.profile(sketches_by_c[c], databases_by_c[c], threads, raw_tsv,
                         estimate_unknown=estimate_unknown)
            raw_tsvs.append(raw_tsv)
        merged = self._merge_raw_profiles(raw_tsvs, os.path.join(working_directory, 'sylph_profile_merged.tsv'))
        if estimate_unknown:
            self._check_unknown_corrected(merged)
        return self.annotate(merged, metapackage, output_annotated_tsv)

    def run_from_sketch(self, sketch_path, metapackage, threads, output_annotated_tsv,
                        working_directory, estimate_unknown=False):
        '''Profile previously-saved sketch(es) against each of the metapackage's
        sylph databases (matching sketch to database by -c), merge, and annotate.'''
        databases = metapackage.sylph_databases()
        if len(databases) == 0:
            raise Exception("Metapackage bundles no {} databases".format(self.BINARY))
        databases_by_c = self._group_databases_by_c(databases)
        raw_tsvs = []
        for i, c in enumerate(sorted(databases_by_c)):
            sketches = self._sketches_for_c(sketch_path, c)
            raw_tsv = os.path.join(working_directory, 'sylph_profile_{}.tsv'.format(i))
            self.profile(sketches, databases_by_c[c], threads, raw_tsv,
                         estimate_unknown=estimate_unknown)
            raw_tsvs.append(raw_tsv)
        merged = self._merge_raw_profiles(raw_tsvs, os.path.join(working_directory, 'sylph_profile_merged.tsv'))
        if estimate_unknown:
            self._check_unknown_corrected(merged)
        return self.annotate(merged, metapackage, output_annotated_tsv)

    def _group_databases_by_c(self, databases):
        '''Group (db_path, c) tuples into {c: [db_path, ...]} so databases sharing
        a -c can be profiled together.'''
        databases_by_c = {}
        for db, c in databases:
            databases_by_c.setdefault(c, []).append(db)
        return databases_by_c

    def _save_sketches_by_c(self, sketches_by_c, sketch_output):
        '''Save sketches into sketch_output/c<C>/ subdirectories, so they can be
        matched back to each database's -c when reused by renew.'''
        os.makedirs(sketch_output, exist_ok=True)
        total = 0
        for c, sketches in sketches_by_c.items():
            subdir = os.path.join(sketch_output, 'c{}'.format(c))
            os.makedirs(subdir, exist_ok=True)
            for s in sketches:
                shutil.copy(s, os.path.join(subdir, os.path.basename(s)))
                total += 1
        logging.info("Saved {} {} sketch(es) to {}/".format(total, self.BINARY, sketch_output))

    def _sketches_for_c(self, sketch_path, c):
        '''Locate the saved sketch(es) made at -c. Accepts a directory written by
        _save_sketches_by_c (c<C>/ subdirs), a flat directory of sketches, or a
        single sketch file.'''
        if os.path.isdir(sketch_path):
            subdir = os.path.join(sketch_path, 'c{}'.format(c))
            search_dir = subdir if os.path.isdir(subdir) else sketch_path
            sketches = sorted(glob.glob(os.path.join(search_dir, '*' + self.SKETCH_EXTENSION)))
        else:
            sketches = [sketch_path]
        if len(sketches) == 0:
            raise Exception("No {} sketch ({}) for c={} found at {}".format(
                self.BINARY, self.SKETCH_EXTENSION, c, sketch_path))
        return sketches

    def _merge_raw_profiles(self, raw_tsvs, output_tsv):
        '''Concatenate raw sylph profile TSVs (one header, then all data rows).'''
        if len(raw_tsvs) == 1:
            return raw_tsvs[0]
        with open(output_tsv, 'w') as out:
            wrote_header = False
            for raw in raw_tsvs:
                with open(raw) as f:
                    header = f.readline()
                    if not header:
                        continue
                    if not wrote_header:
                        out.write(header)
                        wrote_header = True
                    shutil.copyfileobj(f, out)
        return output_tsv

    def _check_unknown_corrected(self, raw_profile_tsv):
        '''Fail loudly if -u did not produce the unknown-corrected coverage it was
        asked for. Callers set alpha to 1 on the strength of that correction, so a
        silently uncorrected profile would be read on the wrong scale.'''
        with open(raw_profile_tsv) as f:
            fields = csv.DictReader(f, delimiter='\t').fieldnames or []
        if 'True_cov' not in fields:
            raise Exception(
                "{} was run with -u (--estimate-unknown) but its profile reports no True_cov "
                "column (found {}). Its coverages are then not on SingleM's scale and cannot "
                "be used with alpha=1.".format(self.BINARY, fields))

    def _extract_accession(self, genome_file):
        match = _GENOME_ACCESSION_REGEX.search(genome_file)
        return match.group(1) if match else None
