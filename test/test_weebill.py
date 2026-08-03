#!/usr/bin/env python3

import unittest
import os.path
import sys
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path
from singlem import weebill
from singlem.sylph import SylphProfiler
from singlem.weebill import WeebillProfiler, is_two_stage_database, read_profiler_for_metapackage


class _StubMetapackage:
    '''Stands in for a Metapackage's sylph_databases and
    genome_accession_to_taxonomy.'''
    def __init__(self, databases, accession_to_taxonomy=None):
        self._databases = databases
        self._map = accession_to_taxonomy or {}

    def sylph_databases(self):
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

    def setUp(self):
        self._original_run = weebill.extern.run

    def tearDown(self):
        weebill.extern.run = self._original_run

    def test_is_two_stage_database(self):
        self.assertTrue(is_two_stage_database('/db/gtdb.syl2db'))
        self.assertFalse(is_two_stage_database('/db/gtdb.syldb'))

    def test_read_profiler_for_metapackage(self):
        # No bundled database at all, so nothing to run.
        self.assertIsNone(read_profiler_for_metapackage(_StubMetapackage([])))
        self.assertIsInstance(
            read_profiler_for_metapackage(_StubMetapackage([('/db/gtdb.syldb', 200)])), SylphProfiler)
        two_stage = read_profiler_for_metapackage(_StubMetapackage([('/db/gtdb.syl2db', 100)]))
        self.assertIsInstance(two_stage, WeebillProfiler)
        self.assertEqual('weebill', two_stage.BINARY)

    def test_read_profiler_for_metapackage_mixed_croaks(self):
        # Neither binary reads both formats, so a mix cannot be profiled.
        with self.assertRaises(Exception):
            read_profiler_for_metapackage(
                _StubMetapackage([('/db/gtdb.syl2db', 100), ('/db/other.syldb', 100)]))

    def test_run_from_reads_profiles_reads_directly(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        metapackage = _StubMetapackage(
            [('/db/gtdb.syl2db', 100)], {'GCF_000744455.1': self.TAXONOMY})
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            out = os.path.join(d, 'annotated.tsv')
            WeebillProfiler().run_from_reads(
                ['r.1.fq.gz'], ['r.2.fq.gz'], metapackage, 4, out, d)
            lines = open(out).read().strip().split('\n')

        # No sketches were wanted on disk, so profiling reads the fastqs directly:
        # one command, not a sketch followed by a profile.
        self.assertEqual(1, len(runner.commands))
        cmd = runner.commands[0]
        self.assertIn('weebill profile --two-stage', cmd)
        self.assertIn(' -u ', cmd + ' ')  # -u, so coverages are on SingleM's scale
        self.assertIn(' -c 100 ', cmd)  # the -c the database was built at
        self.assertIn(' -t 4 ', cmd)
        self.assertIn(' -1 r.1.fq.gz -2 r.2.fq.gz ', cmd)
        self.assertTrue(cmd.endswith('/db/gtdb.syl2db'))

        self.assertEqual('Sample_file\ttaxonomy\tTrue_cov', lines[0])
        self.assertEqual('mock.fq\t{}\t9.8'.format(self.TAXONOMY), lines[1])

    def test_run_from_reads_single_ended(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        metapackage = _StubMetapackage(
            [('/db/gtdb.syl2db', 100)], {'GCF_000744455.1': self.TAXONOMY})
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            WeebillProfiler().run_from_reads(
                ['r.fq.gz'], None, metapackage, 1, os.path.join(d, 'annotated.tsv'), d)
        self.assertIn(' -r r.fq.gz ', runner.commands[0])

    def test_run_from_reads_saving_sketches_sketches_first(self):
        runner = _RecordingRunner()
        weebill.extern.run = runner
        # The sketch command is expected to leave sketches behind; extern.run is
        # stubbed, so write one where weebill would have.
        metapackage = _StubMetapackage(
            [('/db/gtdb.syl2db', 100)], {'GCF_000744455.1': self.TAXONOMY})
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            sketch_dir = os.path.join(d, 'sketch_c100')
            os.makedirs(sketch_dir)
            open(os.path.join(sketch_dir, 'r.sylspc'), 'w').close()
            sketch_output = os.path.join(d, 'saved_sketches')
            WeebillProfiler().run_from_reads(
                ['r.1.fq.gz'], ['r.2.fq.gz'], metapackage, 1,
                os.path.join(d, 'annotated.tsv'), d, sketch_output=sketch_output)
            # Sketches asked for on disk are saved under their -c, so renew can
            # match each back to the database it was made for.
            self.assertTrue(os.path.exists(os.path.join(sketch_output, 'c100', 'r.sylspc')))

        self.assertEqual(2, len(runner.commands))
        self.assertIn('weebill sketch -c 100', runner.commands[0])
        self.assertIn('--compressed-database', runner.commands[0])
        self.assertIn('weebill profile --two-stage', runner.commands[1])
        self.assertIn('.sylspc', runner.commands[1])
        # -c is baked into the sketch, so it is not repeated at profile time.
        self.assertNotIn(' -c ', runner.commands[1])

    def test_run_from_reads_croaks_without_unknown_correction(self):
        # alpha is taken as 1 on the strength of -u, so a profile that came back
        # without the correction must not be passed off as one that has it.
        weebill.extern.run = _RecordingRunner(coverage_column='Eff_cov')
        metapackage = _StubMetapackage([('/db/gtdb.syl2db', 100)])
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            with self.assertRaises(Exception) as context:
                WeebillProfiler().run_from_reads(
                    ['r.1.fq.gz'], None, metapackage, 1, os.path.join(d, 'annotated.tsv'), d)
            self.assertIn('True_cov', str(context.exception))

    def test_run_from_reads_without_database_croaks(self):
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            with self.assertRaises(Exception):
                WeebillProfiler().run_from_reads(
                    ['r.1.fq.gz'], None, _StubMetapackage([]), 1,
                    os.path.join(d, 'annotated.tsv'), d)

    def test_sketches_for_c_uses_compressed_extension(self):
        p = WeebillProfiler()
        with tempfile.TemporaryDirectory(prefix='singlem-weebill-test') as d:
            os.makedirs(os.path.join(d, 'c100'))
            open(os.path.join(d, 'c100', 'x.sylspc'), 'w').close()
            self.assertEqual(['x.sylspc'], [os.path.basename(s) for s in p._sketches_for_c(d, 100)])


if __name__ == "__main__":
    unittest.main()
