#!/usr/bin/env python3

import unittest
import os.path
import sys
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path
from singlem import weebill
from singlem.weebill import WeebillProfiler


class _StubMetapackage:
    '''Stands in for a Metapackage's weebill_databases and
    genome_accession_to_taxonomy.'''
    def __init__(self, databases, accession_to_taxonomy=None):
        self._databases = databases
        self._map = accession_to_taxonomy or {}

    def weebill_databases(self):
        return list(self._databases)

    def genome_accession_to_taxonomy(self, accessions=None):
        if accessions is None:
            return dict(self._map)
        return {a: t for a, t in self._map.items() if a in accessions}


class _RecordingRunner:
    '''Replaces extern.run so the commands built are inspectable, and writes the
    output TSV each command is asked for so the caller can carry on.'''
    def __init__(self, coverage_column='True_cov'):
        self.commands = []
        self.coverage_column = coverage_column

    def __call__(self, cmd):
        self.commands.append(cmd)
        if ' -o ' in cmd:
            output = cmd.split(' -o ')[1].split(' ')[0]
            with open(output, 'w') as f:
                f.write("Sample_file\tGenome_file\t{}\n".format(self.coverage_column))
                f.write("mock.fq\t/db/GCF_000744455.1_genomic.fna.gz\t9.8\n")


class Tests(unittest.TestCase):
    maxDiff = None

    TAXONOMY = 'd__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobacterium_B;s__Methanobacterium_B sp000744455'
    DATABASES = [('/db/gtdb.syl2db', 100)]

    def setUp(self):
        self._original_run = weebill.extern.run

    def tearDown(self):
        weebill.extern.run = self._original_run

    def _metapackage(self):
        return _StubMetapackage(self.DATABASES, {'GCF_000744455.1': self.TAXONOMY})

    def test_extract_accession(self):
        p = WeebillProfiler()
        self.assertEqual('GCF_000744455.1', p._extract_accession(
            'gtdb_genomes_reps_r232/database/GCF/000/744/455/GCF_000744455.1_genomic.fna.gz'))
        self.assertEqual('GCA_000309865.1', p._extract_accession('GCA_000309865.1_genomic.fna'))
        self.assertIsNone(p._extract_accession('not_an_accession.fna'))

    def test_run_from_reads_profiles_reads_directly(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            out = os.path.join(d, 'annotated.tsv')
            WeebillProfiler().run_from_reads(
                ['r.1.fq.gz'], ['r.2.fq.gz'], self._metapackage(), 4, out, d)
            with open(out) as f:
                lines = f.read().strip().split('\n')

        # No sketches were wanted on disk, so profiling reads the fastqs directly:
        # one command, not a sketch followed by a profile.
        self.assertEqual(1, len(runner.commands))
        cmd = runner.commands[0]
        self.assertIn('weebill profile --two-stage', cmd)
        self.assertIn(' -u ', cmd)  # -u, so coverages are on SingleM's scale
        self.assertIn(' -c 100 ', cmd)  # the -c the database was built at
        self.assertIn(' -t 4 ', cmd)
        self.assertTrue(cmd.endswith(' -1 r.1.fq.gz -2 r.2.fq.gz'))
        self.assertIn('/db/gtdb.syl2db', cmd)

        # The unknown-corrected column is carried through annotation, since it is
        # the only record that -u was in effect.
        self.assertEqual('Sample_file\ttaxonomy\tTrue_cov', lines[0])
        self.assertEqual('mock.fq\t{}\t9.8'.format(self.TAXONOMY), lines[1])

    def test_run_from_reads_single_ended(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            WeebillProfiler().run_from_reads(
                ['r.fq.gz'], None, self._metapackage(), 1, os.path.join(d, 'annotated.tsv'), d)
        self.assertTrue(runner.commands[0].endswith(' -r r.fq.gz'))

    def test_run_from_reads_saving_sketches_sketches_first(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            # extern.run is stubbed, so write a sketch where weebill would have.
            sketch_directory = os.path.join(d, 'sketch')
            os.makedirs(sketch_directory)
            open(os.path.join(sketch_directory, 'r.sylspc'), 'w').close()
            sketch_output = os.path.join(d, 'saved_sketches')
            WeebillProfiler().run_from_reads(
                ['r.1.fq.gz'], ['r.2.fq.gz'], self._metapackage(), 1,
                os.path.join(d, 'annotated.tsv'), d, sketch_output=sketch_output)
            self.assertTrue(os.path.exists(os.path.join(sketch_output, 'r.sylspc')))

        self.assertEqual(2, len(runner.commands))
        self.assertIn('weebill sketch -c 100', runner.commands[0])
        self.assertIn('--compressed-database', runner.commands[0])
        self.assertIn('weebill profile --two-stage', runner.commands[1])
        self.assertIn('.sylspc', runner.commands[1])
        # -c is baked into the sketch, so it is not repeated at profile time.
        self.assertNotIn(' -c ', runner.commands[1])

    def test_run_from_sketch(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            saved = os.path.join(d, 'saved_sketches')
            os.makedirs(saved)
            open(os.path.join(saved, 'r.sylspc'), 'w').close()
            WeebillProfiler().run_from_sketch(
                saved, self._metapackage(), 2, os.path.join(d, 'annotated.tsv'), d)
        cmd = runner.commands[0]
        self.assertIn('weebill profile --two-stage -u -t 2', cmd)
        self.assertIn('/db/gtdb.syl2db', cmd)
        self.assertTrue(cmd.endswith('r.sylspc'))

    def test_run_from_sketch_without_sketches_croaks(self):
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            empty = os.path.join(d, 'empty')
            os.makedirs(empty)
            with self.assertRaises(Exception):
                WeebillProfiler().run_from_sketch(
                    empty, self._metapackage(), 1, os.path.join(d, 'annotated.tsv'), d)

    def test_run_from_reads_croaks_without_unknown_correction(self):
        # alpha is taken as 1 on the strength of -u, so a profile that came back
        # without the correction must not be passed off as one that has it.
        weebill.extern.run = _RecordingRunner(coverage_column='Eff_cov')
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            with self.assertRaises(Exception) as context:
                WeebillProfiler().run_from_reads(
                    ['r.1.fq.gz'], None, self._metapackage(), 1,
                    os.path.join(d, 'annotated.tsv'), d)
            self.assertIn('True_cov', str(context.exception))

    def test_run_from_reads_without_database_croaks(self):
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            with self.assertRaises(Exception):
                WeebillProfiler().run_from_reads(
                    ['r.1.fq.gz'], None, _StubMetapackage([]), 1,
                    os.path.join(d, 'annotated.tsv'), d)

    def test_annotate_drops_genomes_the_metapackage_does_not_know(self):
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            raw = os.path.join(d, 'raw.tsv')
            with open(raw, 'w') as f:
                f.write("Sample_file\tGenome_file\tTaxonomic_abundance\tTrue_cov\n")
                f.write("mock.fq\t/db/GCF/GCF_000744455.1_genomic.fna.gz\t95.0\t9.8\n")
                f.write("mock.fq\t/db/GCF/GCF_999999999.9_genomic.fna.gz\t5.0\t1.0\n")
            out = os.path.join(d, 'annotated.tsv')
            WeebillProfiler().annotate(raw, self._metapackage(), out)
            with open(out) as f:
                lines = f.read().strip().split('\n')
        self.assertEqual(2, len(lines))  # header + the one matched row
        self.assertEqual('mock.fq\t{}\t9.8'.format(self.TAXONOMY), lines[1])


if __name__ == "__main__":
    unittest.main()
