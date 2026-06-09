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
    --sylph-profile path consumes (columns Sample_file, taxonomy, Eff_cov).

    Requires the sylph binary on the PATH.'''

    def sketch_reads(self, forward_reads, reverse_reads, c, threads, output_directory):
        '''Sketch reads into output_directory; return the sorted list of .sylsp
        paths produced (one per sample). reverse_reads may be None for single-end.'''
        if not forward_reads:
            raise Exception("No reads provided to sylph sketch")
        cmd = "sylph sketch -c {} -t {} -d {}".format(c, threads, output_directory)
        if reverse_reads:
            cmd += " -1 {} -2 {}".format(' '.join(forward_reads), ' '.join(reverse_reads))
        else:
            cmd += " -r {}".format(' '.join(forward_reads))
        logging.info("Sketching reads with sylph ..")
        extern.run(cmd)
        sylsps = sorted(glob.glob(os.path.join(output_directory, '*.sylsp')))
        if len(sylsps) == 0:
            raise Exception("sylph sketch produced no .sylsp file")
        return sylsps

    def profile(self, sylsp_paths, sylph_dbs, threads, output_tsv):
        '''Profile sample sketches against one or more sylph databases (all of
        which must share the -c baked into the sketches), writing the raw sylph
        TSV to output_tsv. Passing the databases together lets sylph reassign
        shared k-mers across them.'''
        logging.info("Profiling {} sample sketch(es) against {} sylph database(s) ..".format(
            len(sylsp_paths), len(sylph_dbs)))
        extern.run("sylph profile {} {} -t {} -o {}".format(
            ' '.join(sylph_dbs), ' '.join(sylsp_paths), threads, output_tsv))
        return output_tsv

    def annotate(self, raw_profile_tsv, metapackage, output_tsv):
        '''Annotate a raw sylph profile with GTDB taxonomy from the metapackage,
        writing Sample_file / taxonomy / Eff_cov to output_tsv.'''
        rows = []
        needed = set()
        with open(raw_profile_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames is None or 'Genome_file' not in reader.fieldnames or 'Eff_cov' not in reader.fieldnames:
                raise Exception("Unexpected sylph profile format: {}".format(reader.fieldnames))
            for row in reader:
                accession = self._extract_accession(row['Genome_file'])
                rows.append((row.get('Sample_file', ''), accession, row['Eff_cov'], row['Genome_file']))
                if accession is not None:
                    needed.add(accession)

        accession_to_taxonomy = metapackage.genome_accession_to_taxonomy(needed)

        num_written = 0
        num_skipped = 0
        with open(output_tsv, 'w') as out:
            out.write('Sample_file\ttaxonomy\tEff_cov\n')
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
                       output_annotated_tsv, working_directory, sketch_output=None):
        '''Sketch reads (once per distinct -c across the metapackage's sylph
        databases), optionally save the sketches, profile each database against
        the sketch made at its -c, merge, and annotate. Returns output_annotated_tsv.'''
        databases = metapackage.sylph_databases()
        if len(databases) == 0:
            raise Exception("Metapackage bundles no sylph databases")

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
            self.profile(sketches_by_c[c], databases_by_c[c], threads, raw_tsv)
            raw_tsvs.append(raw_tsv)
        merged = self._merge_raw_profiles(raw_tsvs, os.path.join(working_directory, 'sylph_profile_merged.tsv'))
        return self.annotate(merged, metapackage, output_annotated_tsv)

    def run_from_sketch(self, sketch_path, metapackage, threads, output_annotated_tsv, working_directory):
        '''Profile previously-saved sketch(es) against each of the metapackage's
        sylph databases (matching sketch to database by -c), merge, and annotate.'''
        databases = metapackage.sylph_databases()
        if len(databases) == 0:
            raise Exception("Metapackage bundles no sylph databases")
        databases_by_c = self._group_databases_by_c(databases)
        raw_tsvs = []
        for i, c in enumerate(sorted(databases_by_c)):
            sylsps = self._sketches_for_c(sketch_path, c)
            raw_tsv = os.path.join(working_directory, 'sylph_profile_{}.tsv'.format(i))
            self.profile(sylsps, databases_by_c[c], threads, raw_tsv)
            raw_tsvs.append(raw_tsv)
        merged = self._merge_raw_profiles(raw_tsvs, os.path.join(working_directory, 'sylph_profile_merged.tsv'))
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
        for c, sylsps in sketches_by_c.items():
            subdir = os.path.join(sketch_output, 'c{}'.format(c))
            os.makedirs(subdir, exist_ok=True)
            for s in sylsps:
                shutil.copy(s, os.path.join(subdir, os.path.basename(s)))
                total += 1
        logging.info("Saved {} sylph sketch(es) to {}/".format(total, sketch_output))

    def _sketches_for_c(self, sketch_path, c):
        '''Locate the saved sketch(es) made at -c. Accepts a directory written by
        _save_sketches_by_c (c<C>/ subdirs), a flat directory of .sylsp, or a
        single .sylsp file.'''
        if os.path.isdir(sketch_path):
            subdir = os.path.join(sketch_path, 'c{}'.format(c))
            search_dir = subdir if os.path.isdir(subdir) else sketch_path
            sylsps = sorted(glob.glob(os.path.join(search_dir, '*.sylsp')))
        else:
            sylsps = [sketch_path]
        if len(sylsps) == 0:
            raise Exception("No sylph sketch (.sylsp) for c={} found at {}".format(c, sketch_path))
        return sylsps

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

    def _extract_accession(self, genome_file):
        match = _GENOME_ACCESSION_REGEX.search(genome_file)
        return match.group(1) if match else None
