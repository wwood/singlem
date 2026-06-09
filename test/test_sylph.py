#!/usr/bin/env python3

import unittest
import os.path
import sys
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path
from singlem.sylph import SylphProfiler


class _StubMetapackage:
    '''Stands in for a Metapackage's genome_accession_to_taxonomy.'''
    def __init__(self, accession_to_taxonomy):
        self._map = accession_to_taxonomy

    def genome_accession_to_taxonomy(self, accessions=None):
        if accessions is None:
            return dict(self._map)
        return {a: t for a, t in self._map.items() if a in accessions}


class Tests(unittest.TestCase):
    maxDiff = None

    def test_extract_accession(self):
        p = SylphProfiler()
        self.assertEqual('GCF_000744455.1', p._extract_accession(
            'gtdb_genomes_reps_r232/database/GCF/000/744/455/GCF_000744455.1_genomic.fna.gz'))
        self.assertEqual('GCA_000309865.1', p._extract_accession('GCA_000309865.1_genomic.fna'))
        self.assertIsNone(p._extract_accession('not_an_accession.fna'))

    def test_annotate(self):
        tax = 'd__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobacterium_B;s__Methanobacterium_B sp000744455'
        with tempfile.TemporaryDirectory(prefix='singlem-sylph-test') as d:
            raw = os.path.join(d, 'raw.tsv')
            with open(raw, 'w') as f:
                f.write("Sample_file\tGenome_file\tTaxonomic_abundance\tEff_cov\n")
                f.write("mock.fq\t/db/GCF/GCF_000744455.1_genomic.fna.gz\t95.0\t9.8\n")
                # An accession the metapackage doesn't know -> dropped.
                f.write("mock.fq\t/db/GCF/GCF_999999999.9_genomic.fna.gz\t5.0\t1.0\n")
            stub = _StubMetapackage({'GCF_000744455.1': tax})
            out = os.path.join(d, 'annotated.tsv')
            SylphProfiler().annotate(raw, stub, out)
            with open(out) as f:
                lines = f.read().strip().split('\n')
        self.assertEqual('Sample_file\ttaxonomy\tEff_cov', lines[0])
        self.assertEqual(2, len(lines))  # header + one matched row
        self.assertEqual('mock.fq\t{}\t9.8'.format(tax), lines[1])


    def test_merge_raw_profiles(self):
        with tempfile.TemporaryDirectory(prefix='singlem-sylph-test') as d:
            a = os.path.join(d, 'a.tsv')
            b = os.path.join(d, 'b.tsv')
            with open(a, 'w') as f:
                f.write("Sample_file\tGenome_file\tEff_cov\ns\tg1\t5\n")
            with open(b, 'w') as f:
                f.write("Sample_file\tGenome_file\tEff_cov\ns\tg2\t3\ns\tg3\t1\n")
            out = os.path.join(d, 'merged.tsv')
            SylphProfiler()._merge_raw_profiles([a, b], out)
            lines = open(out).read().strip().split('\n')
        self.assertEqual('Sample_file\tGenome_file\tEff_cov', lines[0])
        self.assertEqual(4, len(lines))  # one header + three data rows
        self.assertEqual(['s\tg1\t5', 's\tg2\t3', 's\tg3\t1'], lines[1:])

    def test_merge_raw_profiles_single_passthrough(self):
        with tempfile.TemporaryDirectory(prefix='singlem-sylph-test') as d:
            a = os.path.join(d, 'a.tsv')
            open(a, 'w').write("h\nx\n")
            # A single input is returned as-is (no merge file needed).
            self.assertEqual(a, SylphProfiler()._merge_raw_profiles([a], os.path.join(d, 'm.tsv')))

    def test_sketches_for_c(self):
        p = SylphProfiler()
        with tempfile.TemporaryDirectory(prefix='singlem-sylph-test') as d:
            # Directory written by _save_sketches_by_c: c200/ and c100/ subdirs.
            os.makedirs(os.path.join(d, 'c200'))
            os.makedirs(os.path.join(d, 'c100'))
            open(os.path.join(d, 'c200', 'x.sylsp'), 'w').close()
            open(os.path.join(d, 'c100', 'y.sylsp'), 'w').close()
            self.assertEqual(['x.sylsp'], [os.path.basename(s) for s in p._sketches_for_c(d, 200)])
            self.assertEqual(['y.sylsp'], [os.path.basename(s) for s in p._sketches_for_c(d, 100)])
            # A single .sylsp file is used directly.
            single = os.path.join(d, 'c200', 'x.sylsp')
            self.assertEqual([single], p._sketches_for_c(single, 200))


if __name__ == "__main__":
    unittest.main()
