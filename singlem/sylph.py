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

    def profile(self, sylsp_paths, sylph_db, threads, output_tsv):
        '''Profile sample sketches against the sylph database, writing the raw
        sylph TSV to output_tsv. The -c is baked into the sketches.'''
        logging.info("Profiling {} sample sketch(es) against the sylph database ..".format(len(sylsp_paths)))
        extern.run("sylph profile {} {} -t {} -o {}".format(
            sylph_db, ' '.join(sylsp_paths), threads, output_tsv))
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
        '''Sketch reads, optionally save the sketch(es), profile against the
        metapackage's sylph DB, and annotate. Returns output_annotated_tsv.'''
        sketch_dir = os.path.join(working_directory, 'sylph_sketch')
        os.makedirs(sketch_dir, exist_ok=True)
        sylsps = self.sketch_reads(
            forward_reads, reverse_reads, metapackage.sylph_c(), threads, sketch_dir)
        if sketch_output is not None:
            self._save_sketches(sylsps, sketch_output)
        raw_tsv = os.path.join(working_directory, 'sylph_profile.tsv')
        self.profile(sylsps, metapackage.sylph_db_path(), threads, raw_tsv)
        return self.annotate(raw_tsv, metapackage, output_annotated_tsv)

    def run_from_sketch(self, sketch_path, metapackage, threads, output_annotated_tsv, working_directory):
        '''Profile a previously-saved sketch (a .sylsp file or a directory of
        them) against the metapackage's sylph DB, and annotate.'''
        if os.path.isdir(sketch_path):
            sylsps = sorted(glob.glob(os.path.join(sketch_path, '*.sylsp')))
        else:
            sylsps = [sketch_path]
        if len(sylsps) == 0:
            raise Exception("No .sylsp sketch found at {}".format(sketch_path))
        raw_tsv = os.path.join(working_directory, 'sylph_profile.tsv')
        self.profile(sylsps, metapackage.sylph_db_path(), threads, raw_tsv)
        return self.annotate(raw_tsv, metapackage, output_annotated_tsv)

    def _save_sketches(self, sylsps, sketch_output):
        '''Save a single sketch to the given path, or multiple into a directory.'''
        if len(sylsps) == 1:
            shutil.copy(sylsps[0], sketch_output)
            logging.info("Saved sylph sketch to {}".format(sketch_output))
        else:
            os.makedirs(sketch_output, exist_ok=True)
            for s in sylsps:
                shutil.copy(s, os.path.join(sketch_output, os.path.basename(s)))
            logging.info("Saved {} sylph sketches to {}/".format(len(sylsps), sketch_output))

    def _extract_accession(self, genome_file):
        match = _GENOME_ACCESSION_REGEX.search(genome_file)
        return match.group(1) if match else None
