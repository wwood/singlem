#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft, Tim Lamberton.
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import os.path
import tempfile
import extern
import sys
import json
import re

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.pipe import SearchPipe
from singlem.sequence_classes import SeqReader

TEST_ANNOY = False
try:
    import annoy
    TEST_ANNOY = True
    print("annoy found, running relevant tests", file=sys.stderr)
except ImportError:
    print("WARNING: annoy not found, skipping relevant tests", file=sys.stderr)
    pass

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    headers_with_extras = headers + str.split('read_names nucleotides_aligned taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
    maxDiff = None
    two_packages = '%s %s' % (
        os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
        os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def assertEqualOtuTable(self, expected_array_or_string, observed_string, no_assign_taxonomy=False):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])

        r = re.compile(r'  +')
        if isinstance(expected_array_or_string, str):
            expected_array = list(expected_array_or_string.split("\n"))
            expected_array = [r.sub("\t", line) for line in expected_array]
            expected_array = [line.split("\t") for line in expected_array]
            if no_assign_taxonomy:
                expected_array2 = [expected_array[0]]
                for line in expected_array[1:]:
                    expected_array2.append(line+[''])
                expected_array = expected_array2
        else:
            expected_array = expected_array_or_string

        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))
    
    def test_fast_protein_package(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --no-diamond-prefilter --assignment-method diamond --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --assignment-method diamond --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_diamond_package_assignment(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --assignment-method diamond --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_diamond_package_assignment_multithreaded(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            # By specifying --threads 2 multiprocessing is used in package assignment, even if there is only 1 package
            cmd = "%s pipe --threads 2 --sequences %s --assignment-method diamond --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_diamond_package_assignment_paired_both_hit(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs)
                n2.flush()

                cmd = "%s pipe --forward %s --reverse %s --assignment-method diamond --otu-table /dev/stdout --singlem-packages %s" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_diamond_package_assignment_paired_one_hits(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write('''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA''')
                n2.flush()

                cmd = "%s pipe --assignment-method diamond --forward %s --reverse %s --otu-table /dev/stdout --singlem-packages %s" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))
                cmd = "%s pipe --assignment-method diamond --reverse %s --forward %s --otu-table /dev/stdout --singlem-packages %s" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n2.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --assignment-method diamond --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_prefilter_with_diamond_assignment(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s --assignment-method diamond" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter_with_diamond_assignment_metapackage(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --metapackage %s --assignment-method diamond" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg.smpkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_prefilter_with_diamond_assignment_metapackage_condensed_outputs(self):
        expected = \
            "sample\tcoverage\ttaxonomy\n" + \
            'SAMPLE_ID	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]\n'
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --metapackage %s --assignment-method diamond --taxonomic-profile /dev/stdout --taxonomic-profile-krona /tmp/blah.html" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg.smpkg'))
            print(cmd)
            self.assertEqual(
                expected, extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),'SAMPLE_ID'))

    def test_fast_protein_package_prefilter_with_diamond_assignment_2_samples(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='sample1.fa') as n:
            n.write(inseqs)
            n.flush()

            with tempfile.NamedTemporaryFile(mode='w',suffix='sample2.fa') as n2:
                n2.write(inseqs)
                n2.flush()

                cmd = "%s pipe --sequences %s %s --otu-table /dev/stdout --singlem-packages %s --assignment-method diamond" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),'').replace(os.path.basename(n2.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter_with_diamond_assignment_paired(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --forward %s --reverse %s --otu-table /dev/stdout --singlem-packages %s --assignment-method diamond" % (
                path_to_script, n.name, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_prefilter_with_diamond_example_assignment(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	2524614704',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --forward %s --reverse %s --otu-table /dev/stdout --singlem-packages %s --assignment-method diamond_example" % (
                path_to_script, n.name, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))
    
    def test_minimal(self):
        expected = [
            self.headers,
            ['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --sequences %s/1_pipe/minimal.fa --assignment-method diamond --otu-table /dev/stdout --threads 4 --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg" % (path_to_script, path_to_data, path_to_data)
        
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))
        
    def test_minimal_no_prefilter(self):
        expected = [
            self.headers,
            ['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --sequences %s/1_pipe/minimal.fa --assignment-method diamond --otu-table /dev/stdout --threads 4 --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg --no-diamond-prefilter" % (
            path_to_script,
            path_to_data,
            path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_insert(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','2','4.95','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --sequences %s/1_pipe/insert.fna --assignment-method diamond --no-diamond-prefilter --otu-table /dev/stdout --threads 4 --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg" % (
            path_to_script,
            path_to_data,
            path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))
    
    def test_insert_prefilter(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','2','4.95','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --sequences %s/1_pipe/insert.fna --assignment-method diamond --otu-table /dev/stdout --threads 4 --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg" % (
            path_to_script,
            path_to_data,
            path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_print_insert(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes'],
                    ['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTtttCAAGCAGGTGTG','1','2.51','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --debug --sequences %s/1_pipe/insert.fna --otu-table /dev/stdout --assignment-method diamond --threads 4 --include-inserts --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg" % (
            path_to_script,
            path_to_data,
            path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    @unittest.skip("not really supported any more")
    def test_known_tax_table(self):
        expected = [
            self.headers,
            ['4.12.22seqs','small',
             'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
             '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Gracilibacillus; s__Gracilibacillus_lacisalsi'],
            ['4.11.22seqs','small',
             'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
             '2','4.88','Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --quiet --sequences %s/1_pipe/small.fa --no-diamond-prefilter --otu-table /dev/stdout --assignment-method diamond --threads 4 --singlem-packages %s" % (
            path_to_script,
            path_to_data,
            self.two_packages)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

        expected = [self.headers,
                    ['4.12.22seqs','small',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','some1'],
                    ['4.11.22seqs','small',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        with tempfile.NamedTemporaryFile(mode='w',prefix='singlem_test_known') as t:
            t.write('\n'.join(["\t".join(x) for x in expected[:2]]))
            t.flush()

            cmd = "%s --quiet pipe --sequences %s/1_pipe/small.fa --no-diamond-prefilter --otu-table /dev/stdout --threads 4 --known_otu_tables %s --singlem-packages %s"\
                 % (path_to_script,
                    path_to_data,
                    t.name,
                    self.two_packages)
            self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_diamond_assign_taxonomy(self):
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fasta') as f:
            query = "\n".join(['>HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910 1:N:0:TCCTGAGCCTAAGCCT',
                'GTTAAATTACAAATTCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATCATGGGATTCTGTAAAGAGT',''])
            f.write(query)
            f.flush()

            cmd = "%s pipe --debug --sequences %s --otu-table /dev/stdout --assignment_method diamond --threads 4 --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg" % (
                path_to_script,
                f.name,
                path_to_data)

            expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK',os.path.basename(f.name)[:-6],'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes; c__Bacilli_A; o__Thermoactinomycetales; f__Thermoactinomycetaceae']]
            expected = ["\t".join(x) for x in expected]+['']
            observed = extern.run(cmd).split("\n")
            r = re.compile('; c__.*') # Do not test beyond phylum level because updated diamond version change slightly.
            self.assertEqual([r.sub('',e) for e in expected], [r.sub('',e) for e in observed])

    def test_diamond_example_assign_taxonomy(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2513237297']
                    ]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s pipe --debug --sequences %s/1_pipe/minimal.fa --otu-table /dev/stdout --threads 4 --assignment_method diamond_example --metapackage %s/S1.5.ribosomal_protein_L11_rplK.gpkg.spkg.smpkg" % (
            path_to_script,
            path_to_data,
            path_to_data)
        observed = sorted(extern.run(cmd).split("\n"))
        r = re.compile('\t.*?$') # Do not test the exact genome number because updated diamond version change this slightly.
        self.assertEqual([r.sub('',e) for e in exp], [r.sub('',e) for e in observed])

    def test_one_read_two_orfs_two_diamond_hits(self):
        # what a pain the real world is
        seq = '''>HWI-ST1240:128:C1DG3ACXX:7:2204:6599:65352 1:N:0:GTAGAGGATAGATCGC
ACCCACAGCTCGGGGTTGCCCTTGCCCGACCCCATGCGTGTCTCGGCGGGCTTCTGGTGACGGGCTTGTCCGGGAAGACGCGGATCCAGACCTTGCCTCCGCGCTTGACGTGCCGGGTCATCGCGATACGGGCCGCCTCGATCTGACGTGC
'''
        expected = [
            self.headers,
            ['S1.7.ribosomal_protein_L16_L10E_rplP		CGCGTCTTCCCGGACAAGCCCGTCACCAGAAGCCCGCCGAGACACGCATGGGGTCGGGCA	1	1.64	GCA_000949295.1']]
        exp = sorted(["\t".join(x) for x in expected]+[''])
        with tempfile.NamedTemporaryFile(mode='w',prefix='singlem_test',suffix='.fa') as t:
            t.write(seq)
            t.flush()
            cmd = "%s pipe --quiet --sequences %s --otu-table /dev/stdout --threads 4 --assignment_method diamond_example --singlem-packages %s/S1.7.ribosomal_protein_L16_L10E_rplP.gpkg.spkg" % (
                path_to_script,
                t.name,
                path_to_data)
            self.assertEqual(exp,
                             sorted(extern.run(cmd).
                                    replace(
                                        os.path.basename(t.name).replace('.fa',''),
                                        '').
                                    split("\n")))

    @unittest.skip(reason="fails now due to pplacer/jplace versioning issues I think")
    def test_jplace_output(self):
        expected_jpace = {'fields': ['classification',
                                     'distal_length',
                                     'edge_num',
                                     'like_weight_ratio',
                                     'likelihood',
                                     'pendant_length'],
                          'metadata': 'the_metadata',
                          'placements':
                          [{
                              'nm': [['CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                                      2]],
                              'p': [
                                  [
                                      "o__Bacillales",
                                      0.0874346630859,
                                      13,
                                      0.333350512423,
                                      -608.20180926,
                                      6.11351501465e-06
                                  ],
                                  [
                                      "o__Bacillales",
                                      0.0643521435547,
                                      14,
                                      0.333326884837,
                                      -608.201880142,
                                      6.11351501465e-06
                                  ],
                                  [
                                      "p__Firmicutes",
                                      5.97534179688e-06,
                                      15,
                                      0.33332260274,
                                      -608.201892989,
                                      6.11351501465e-06
                                  ]
                              ]}],
                          'tree': 'tree_thanks',
                          'version': 3}

        with tempfile.TemporaryDirectory(prefix='singlem-test') as d:
            cmd = "%s pipe --sequences %s --assignment-method pplacer --otu-table /dev/null --output-jplace %s"\
                  " --singlem-packages %s" % (
                      path_to_script,
                      os.path.join(path_to_data,'1_pipe','jplace_test.fna'),
                      os.path.join(d, "my_jplace"),
                      os.path.join(path_to_data,'4.12.22seqs.spkg'))
            extern.run(cmd)
            jplace_path = os.path.join(d, 'my_jplace_jplace_test_4.12.22seqs.jplace')
            with open(jplace_path) as f:
                j = json.load(f)
            j['tree'] = 'tree_thanks'
            j['metadata'] = 'the_metadata'
            self.assertEqual(expected_jpace, j)

            # Make sure the guppy sing does not croak
            extern.run("guppy sing -o /dev/null '%s'" % jplace_path)

    def test_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	1	2.56	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))
            
    def test_nucleotide_package_prefilter(self):
        """ correct behaviour is to fail, as DIAMOND does not have a blastn capability """
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            try:
                cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s" % (
                    path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
                self.fail() # this is meant to fail
            except:
                pass
                

    def test_revcom_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	1	2.56	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_fwd_and_revcom_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	2	5.13	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_two_nucleotide_packages(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	2	5.13	Root; k__Bacteria; p__Proteobacteria',
            '61_otus.second.v3		TTAGGTAGTTGCTGGGGTAACGTCCCAACAAGCCGATAATCGGTACGGGTTGTGAGAGCA	1	1.66	Root; k__Archaea; p__Euryarchaeota',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
>NS500333:10:H0V2GAGXX:2:13211:8623:16289 1:N:0:GATCAG
ATTAGGTAGTTGCTGGGGTAACGTCCCAACAAGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGATGGATTCTGAGACACGAATCCAGGTCCTACGGGGCGCAGCAGGCGCGAAAACTTTACACTGCGCGAAAGCGCGATA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s %s" % (
                path_to_script,
                n.name,
                os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'),
                os.path.join(path_to_data,'second_packge.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_no_taxonomy(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s --no-assign-taxonomy" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_sample_name_strange_characters(self):
        expected = [self.headers,
                    ['4.12.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.8',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Gracilibacillus; s__Gracilibacillus_lacisalsi'],
                    ['4.11.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.8',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179'],
                    ['4.12.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.9',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Gracilibacillus; s__Gracilibacillus_lacisalsi'],
                    ['4.11.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.9',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']]
        cmd = "%s pipe --quiet --sequences "\
        "%s/1_pipe/contigs.fasta.metabat-bins-_-t20_--superspecific.9.fa "\
        "%s/1_pipe/contigs.fasta.metabat-bins-_-t20_--superspecific.8.fa "\
        "--otu-table /dev/stdout --threads 4 --assignment-method diamond --singlem-packages %s" %(
            path_to_script,
            path_to_data,
            path_to_data,
            self.two_packages)
        self.assertEqualOtuTable(expected, extern.run(cmd))

    def test_archive_otu_groopm_compatibility(self):
        """This tests for API stability, where the API is used by GroopM 2.0"""
        expected = [('contig_1', '4.11.22seqs', 'Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179')]

        inseqs = '''>contig_1 abcd
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --archive-otu-table /dev/stdout --singlem-packages %s --assignment-method diamond" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))

            j = json.loads(extern.run(cmd))
            fields = j['fields']
            data = j['otus']
            self.assertEqual(expected,
                             [(name, row[fields.index('gene')], row[fields.index('taxonomy')]) for row in data for name in row[fields.index('read_names')]]
                            )

    def test_protein_package_non60_length(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTT	1	1.41	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --assignment-method diamond --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.length30.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_nucleotide_package_non60_length(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACAC	1	1.10	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.length10.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_nucleotide_matches_forward_and_reverse(self):
        '''Not the best test since the read does not actually match the alignment
        region, but at least no error should be thrown.'''
        expected = [
            "\t".join(self.headers),
            '']
        inseqs = '''>seq18975201
GTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGAGTCGAGTTGCAGACTCCAATCCGAACTGGGGCCGGTTTTTATGGATTGGCTTCCCCTCGCGGGTTCGCGACCCTTTGTACCGGCCATTGTAACACGTGTGTAGC
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu-table /dev/stdout --singlem-packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_paired_reads_hello_world(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --sequences {} --otu-table /dev/stdout --assignment-method diamond --singlem-packages {} --reverse {}".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_paired_reads_one_read_each(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2	60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT	[\'2524614704\'] [\'2524614704\']	diamond',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --sequences {} --otu-table /dev/stdout --assignment-method diamond --singlem-packages {} --reverse {} --output-extras".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_paired_reads_two_reads_each(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	3	7.32	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2 seq3	60 60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA	[\'2524614704\'] [\'2524614704\'] [\'2524614704\']	diamond',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
>seq3
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
>seq3
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --sequences {} --otu-table /dev/stdout --assignment-method diamond --singlem-packages {} --reverse {} --output-extras".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_paired_reads_one_read_each_diamond(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2	60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT	[\'2524614704\'] [\'2524614704\']	diamond',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --reverse {} --output-extras --assignment_method diamond".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_paired_reads_one_read_each_diamond_example(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	2524614704	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2	60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT	None	diamond_example',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --reverse {} --output-extras --assignment_method diamond_example".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_two_orfs_in_same_read(self):
        # Read finds 2 ORFs, pretty rare for reads (but happens for genomes
        # more frequently)
        expected = [
            "\t".join(self.headers_with_extras),
            'S1.12.ribosomal_protein_S12_S23		CGTGGTGTCTGCACCCGGGTGTACACCACCACCCGAAGA---------AGCCGAACTCGG	1	1.50	Root; d__Bacteria; p__Actinobacteria; c__Acidimicrobiia; o__Microtrichales; f__TK06; g__MedAcidi-G3; s__MedAcidi-G3_sp4	A00178:38:H5NYYDSXX:2:1552:32524:1517	51	False	CGGGATGTAGGCAGTGACCTCCACGCCTGAGGAGAGCCGGACGCGTGCGACCTTGCGCAACGCCGAGTTCGGCTTCTTCGGGTGGTGGTGTACACCCGGGTGCAGACACCACGGCGCTGGGGCGAACCCTTGAGCGCAGGGGTGTTGGTCT	[\'GCA_000817105.1\']	diamond',
            '']
        inseqs = '''>A00178:38:H5NYYDSXX:2:1552:32524:1517 1:N:0:CAACGGA+ATCCGTT
CGGGATGTAGGCAGTGACCTCCACGCCTGAGGAGAGCCGGACGCGTGCGACCTTGCGCAACGCCGAGTTCGGCTTCTTCGGGTGGTGGTGTACACCCGGGTGCAGACACCACGGCGCTGGGGCGAACCCTTGAGCGCAGGGGTGTTGGTCT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --output-extras --assignment-method diamond".format(
                path_to_script,
                n.name,
                os.path.join(path_to_data, 'S1.12.ribosomal_protein_S12_S23.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_split_genes(self):
        # 2 ORFs found causing GraftM to do its "split" thing
        expected = [
            "\t".join(self.headers),
            "S1.2.ribosomal_protein_L3_rplC	aa_orf_split_bug	GTTGACGTGGCGGCCATCACAAAGGGCAAGGGATGGCAGGGCGTCCTGAAGCGGTGGAAC	1	1.05	Root; d__Archaea; p__Crenarchaeota; c__Nitrososphaeria; o__Nitrososphaerales; f__Nitrosopumilaceae; g__Nitrosopumilus; s__Nitrosopumilus_piranensis",
            '']
        cmd = "{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --assignment-method diamond".format(
            path_to_script,
            os.path.join(path_to_data, 'aa_orf_split_bug.fna'),
            os.path.join(path_to_data, 'S1.2.ribosomal_protein_L3_rplC.gpkg.spkg'))
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_no_results(self):
        # ORFs but no hits at all 
        expected = [
            "\t".join(self.headers),
            '']
        cmd = "{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {}".format(
            path_to_script,
            os.path.join(path_to_data, 'random.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'))
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_no_results_diamond_prefilter(self):
        # ORFs but no hits at all 
        expected = [
            "\t".join(self.headers),
            '']
        cmd = "{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {}".format(
            path_to_script,
            os.path.join(path_to_data, 'random.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'))
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_genome(self):
        expected = [
            "\t".join(self.headers),
            'S1.13.ribosomal_protein_S15P_S13e	uap2	AAGGATTTGAGTGCAAAAAGAGGACTCGATTTTACAGAGGCAAAGATAAGAAAACTTGGA	1	1.00	Root; d__Archaea; p__Halobacterota; c__Methanomicrobia; o__Methanomicrobiales; f__Methanocullaceae; g__Methanoculleus',
            '']

        cmd = '{} pipe --genome-fasta-files {}//uap2.fna --singlem-packages {}/S1.13.chainsaw.gpkg.spkg --otu-table /dev/stdout --assignment-method diamond --assignment-method diamond'.format(
            path_to_script,
            path_to_data,
            path_to_data,
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_prefilter_piped_in_input(self):
        with tempfile.NamedTemporaryFile(suffix='.mkfifo.fna') as fifo:
            with tempfile.NamedTemporaryFile(suffix='.mkfifo.fna') as script:
                script.write("mkfifo {}\n".format(fifo.name).encode())
                script.write("cat {}/1_pipe/small.fa > {} &\n".format(path_to_data, fifo.name).encode())
                script.write('{} pipe --forward {} --singlem-packages {}/4.11.22seqs.gpkg.spkg --otu-table /dev/stdout --assignment-method diamond --assignment-method diamond'.format(
                    path_to_script,
                    fifo.name,
                    path_to_data,
                ).encode())
                script.flush()

                expected = [
                    "\t".join(self.headers),
                    '4.11.22seqs	{}	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179'.format(os.path.basename(fifo.name.replace('.fna',''))),
                    '']
                
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run('bash {}'.format(script.name)))

    @unittest.skipIf(not TEST_ANNOY, "annoy not installed")
    def test_annoy_only_assignment_single(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method annoy'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    @unittest.skipIf(not TEST_ANNOY, "annoy not installed")
    def test_annoy_only_assignment_paired(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --forward {} --reverse {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method annoy'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs2.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    @unittest.skipIf(not TEST_ANNOY, "annoy not installed")
    def test_annoy_only_assignment_output_extras_single(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy	read_names	nucleotides_aligned	taxonomy_by_known?	read_unaligned_sequences	equal_best_hit_taxonomies	taxonomy_assignment_method',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482	60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179 Root; part_of_sdb; novel_domain	singlem_query_based']
        cmd = '{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method annoy --output-extras'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs_and_inseqs2.manually_different_species.otu_table.csv.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    @unittest.skipIf(not TEST_ANNOY, "annoy not installed")
    def test_annoy_only_assignment_output_archive_single(self):
        cmd = '{} pipe --sequences {} --archive-otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method annoy'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs_and_inseqs2.manually_different_species.otu_table.csv.sdb'),
        )
        self.assertEqual(json.loads('{"version": 4, "alignment_hmm_sha256s": ["4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9"], "singlem_package_sha256s": ["e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48"], "fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?", "read_unaligned_sequences", "equal_best_hit_taxonomies", "taxonomy_assignment_method"], "otus": [["4.11.22seqs", "4.11.22seqs.gpkg.spkg_inseqs", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 1, 2.4390243902439024, "Root; part_of_sdb", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482"], [60], false, ["ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA"], ["Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179", "Root; part_of_sdb; novel_domain"], "singlem_query_based"]]}'),
            json.loads(extern.run(cmd)))

    def test_smafa_naive_then_diamond_single(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method smafa_naive_then_diamond'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.smafa_naive.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    @unittest.skip('scann not currently working on small DBs')
    def test_scann_then_diamond_assignment_single(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --sequences {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method scann_then_diamond'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_smafa_naive_then_diamond_assignment_paired(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --forward {} --reverse {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method smafa_naive_then_diamond'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs2.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.smafa_naive.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_smafa_naive_then_diamond_assignment_paired_only_read2_hits(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	random	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --forward {} --reverse {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method smafa_naive_then_diamond'.format(
            path_to_script,
            os.path.join(path_to_data, 'random.fna'), # i.e. no hits
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs2.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.smafa_naive.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_smafa_naive_then_diamond_assignment_paired_only_read1_hits(self):
        expected = ['gene	sample	sequence	num_hits	coverage	taxonomy',
        '4.11.22seqs	4.11.22seqs.gpkg.spkg_inseqs2	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; part_of_sdb; Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179']
        cmd = '{} pipe --forward {} --reverse {} --otu-table /dev/stdout --singlem-packages {} --assignment-singlem-db {} --assignment-method smafa_naive_then_diamond'.format(
            path_to_script,
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg_inseqs2.fna'),
            os.path.join(path_to_data, 'random.fna'), # i.e. no hits
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
            os.path.join(path_to_data, '4.11.22seqs.paired.manual.json.v5.smafa_naive.sdb'),
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_exclude_off_target_hits(self):
        without_exclude = "gene    sample  sequence        num_hits        coverage        taxonomy\n" \
                "4.11.22seqs     4.11.22seqs.gpkg.spkg_inseqs    TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA    1       2.44    Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179\n"
        cmd = f"{path_to_script} pipe --forward {path_to_data}/4.11.22seqs.gpkg.spkg_inseqs.fna --assignment-method diamond --singlem-package {path_to_data}/4.11.22seqs.v3_archaea_targetted.gpkg.spkg --otu-table /dev/stdout"
        self.assertEqualOtuTable(without_exclude, extern.run(cmd))

        with_exclude = "gene    sample  sequence        num_hits        coverage        taxonomy\n"
        cmd = f"{path_to_script} pipe --forward {path_to_data}/4.11.22seqs.gpkg.spkg_inseqs.fna --assignment-method diamond --singlem-package {path_to_data}/4.11.22seqs.v3_archaea_targetted.gpkg.spkg --otu-table /dev/stdout --exclude-off-target-hits"
        self.assertEqualOtuTable(with_exclude, extern.run(cmd))

    def test_genome_input_dereplication(self):
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            '4.12.22seqs     GCA_000309865.1_genomic  GATGGCGGTAAAGCCACTCCCGGCCCACCATTAGGTCCAGCAATCGGACCCCTAGGTATC    1       1.00    '
        # ~/git/singlem/bin/singlem pipe --genome-fasta-files genomes/GCA_000309865.1_genomic.fna --singlem-package ../4.12.22seqs.spkg/ --otu-table /dev/stdout --no-assign-taxonomy --min-orf-length 96
        cmd = '{} pipe --translation-table 11 --genome-fasta-files {} --singlem-package {} --otu-table /dev/stdout --no-assign-taxonomy --min-orf-length 96'.format(
            path_to_script,
            os.path.join(path_to_data, 'methanobacteria/genomes/GCA_000309865.1_genomic.fna'),
            os.path.join(path_to_data, '4.12.22seqs.spkg'),
        )
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_genome_multiplexing(self):
        cmd_stub = '{} pipe --metapackage {} --output-extras --otu-table /dev/stdout --assignment-method diamond --genome-fasta-files '.format(
            path_to_script,
            os.path.join(path_to_data, 'archaeal_small.v4.smpkg'),
        )
        g1 = os.path.join(path_to_data, 'methanobacteria/genomes/GCA_000309865.1_genomic.fna')
        g2 = os.path.join(path_to_data, 'methanobacteria/genomes/GCF_000191585.1_genomic.fna')
        cmd_t1 = cmd_stub + g1
        cmd_t2 = cmd_stub + g2
        cmd_both = cmd_stub + g1 + ' ' + g2

        r1 = extern.run(cmd_t1)
        r2 = extern.run(cmd_t2)
        r_both = extern.run(cmd_both)
        
        lines2 = r2.split('\n')
        r2 = '\n'.join(lines2[1:]) # Remove header
        self.assertEqualOtuTable(r1+r2, r_both)

    def test_translation_table4_no_diamond_prefilter(self):
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  tt4_s1.2        ATAAACTTAATAGGTACATCAAAAGGTAAAGGTTTTCAATGAGTTATGAAAAGATTTCAT    1       1.11    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__RF39; f__CAG-1000; g__CAG-1000'
        cmd = f'{path_to_script} pipe --forward {path_to_data}/tt4_s1.2.fna --otu-table /dev/stdout --no-diamond-prefilter --translation-table 4 --threads 32 --singlem-package {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_translation_table4_diamond_prefilter(self):
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  tt4_s1.2        ATAAACTTAATAGGTACATCAAAAGGTAAAGGTTTTCAATGAGTTATGAAAAGATTTCAT    1       1.11    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__RF39; f__CAG-1000; g__CAG-1000'
        cmd = f'{path_to_script} pipe --forward {path_to_data}/tt4_s1.2.fna --otu-table /dev/stdout --translation-table 4 --threads 32 --singlem-package {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra1(self):
        '''
        Run on SRR8653040.sra
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    20       26.15    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTTACAGC    1       1.31    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra_chunk1(self):
        '''
        Run on SRR8653040.sra, which has 424064 reads. This test only runs on the first 200,000 reads.
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    13       17.00    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTTACAGC    1       1.31    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond --read-chunk-number 1 --read-chunk-size 200000'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra_chunk2(self):
        '''
        Run on SRR8653040.sra, which has 424064 reads. This test only runs on the first 200,000 reads.
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    6       7.84   Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond --read-chunk-number 2 --read-chunk-size 200000'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))

    def test_sra_chunk3(self):
        '''
        Run on SRR8653040.sra, which has 424064 reads. This test only runs on the first 200,000 reads.
        '''
        expected = 'gene    sample  sequence        num_hits        coverage        taxonomy\n' \
            'S1.2.ribosomal_protein_L3_rplC  SRR8653040      GTTGATGTTACAGGTACTACGAAAGGTAAAGGATTCCAAGGGGCAATCAAACGTCACGGC    1       1.31    Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_faecalis'
        cmd = f'{path_to_script} pipe --sra {path_to_data}/SRR8653040.sra --otu-table /dev/stdout --singlem-packages {path_to_data}/S1.2.ribosomal_protein_L3_rplC.gpkg.spkg/ --assignment-method diamond --read-chunk-number 3 --read-chunk-size 200000'
        self.assertEqualOtuTable(
            expected,
            extern.run(cmd))




if __name__ == "__main__":
    unittest.main()
