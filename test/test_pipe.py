
#!/usr/bin/env python

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
import subprocess
import os.path
import tempfile
import tempdir
from string import split
import extern
import sys
import json
import re

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.pipe import SearchPipe

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None
    two_packages = '%s %s' % (
        os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
        os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    def test_fast_protein_package(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_minimal(self):
        expected = [
            self.headers,
            ['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --sequences %s/1_pipe/minimal.fa --otu_table /dev/stdout --threads 4" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_insert(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','2','4.95','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --quiet pipe --sequences %s/1_pipe/insert.fna --otu_table /dev/stdout --threads 4" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(subprocess.check_output(cmd, shell=True).split("\n")))

    def test_print_insert(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTtttCAAGCAGGTGTG','1','2.51','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --sequences %s/1_pipe/insert.fna --otu_table /dev/stdout --threads 4 --include_inserts" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_known_tax_table(self):
        expected = [
            self.headers,
            ['4.12.22seqs','small',
             'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
             '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.22seqs','small',
             'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
             '2','4.88','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --quiet pipe --sequences %s/1_pipe/small.fa --otu_table /dev/stdout --threads 4 --singlem_packages %s" % (
            path_to_script,
            path_to_data,
            self.two_packages)
        self.assertEqual(exp, sorted(subprocess.check_output(cmd, shell=True).split("\n")))

        expected = [self.headers,
                    ['4.12.22seqs','small',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','some1'],
                    ['4.11.22seqs','small',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        with tempfile.NamedTemporaryFile(prefix='singlem_test_known') as t:
            t.write('\n'.join(["\t".join(x) for x in expected[:2]]))
            t.flush()

            cmd = "%s --quiet pipe --sequences %s/1_pipe/small.fa --otu_table /dev/stdout --threads 4 --known_otu_tables %s --singlem_packages %s"\
                 % (path_to_script,
                    path_to_data,
                    t.name,
                    self.two_packages)
            self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_diamond_assign_taxonomy(self):
        with tempfile.NamedTemporaryFile(suffix='.fasta') as f:
            query = "\n".join(['>HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910 1:N:0:TCCTGAGCCTAAGCCT',
                'GTTAAATTACAAATTCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATCATGGGATTCTGTAAAGAGT',''])
            f.write(query)
            f.flush()

            cmd = "%s --debug pipe --sequences %s --otu_table /dev/stdout --assignment_method diamond --threads 4" % (path_to_script,
                                                            f.name)

            expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK',os.path.basename(f.name)[:-6],'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes; c__Bacilli_A; o__Thermoactinomycetales; f__Thermoactinomycetaceae']]
            expected = ["\t".join(x) for x in expected]+['']
            observed = extern.run(cmd).split("\n")
            r = re.compile('; g__.*') # Do not test beyond genus level because updated diamond version change slightly.
            self.assertEqual([r.sub('',e) for e in expected], [r.sub('',e) for e in observed])

    def test_diamond_example_assign_taxonomy(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2513237297']
                    ]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --sequences %s/1_pipe/minimal.fa --otu_table /dev/stdout --threads 4 --assignment_method diamond_example" % (path_to_script,
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
        with tempfile.NamedTemporaryFile(prefix='singlem_test',suffix='.fa') as t:
            t.write(seq)
            t.flush()
            cmd = "%s --quiet pipe --sequences %s --otu_table /dev/stdout --threads 4 --assignment_method diamond_example" % (path_to_script,
                                                                                                    t.name)
            self.assertEqual(exp,
                             sorted(extern.run(cmd).
                                    replace(
                                        os.path.basename(t.name).replace('.fa',''),
                                        '').
                                    split("\n")))

    def test_jplace_output(self):
        expected_jpace = {u'fields': [u'classification',
                                      u'distal_length',
                                      u'edge_num',
                                      u'like_weight_ratio',
                                      u'likelihood',
                                      u'pendant_length'],
                          u'metadata': 'the_metadata',
                          u'placements':
                          [{
                           u'nm': [[u'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                                     2]],
                            u'p': [[u'o__Bacillales',
                                    0.0874346630859,
                                    13,
                                    0.333351177694,
                                    -631.301684875,
                                    0.150831104822],
                                   [u'o__Bacillales',
                                    0.0643521435547,
                                    14,
                                    0.333326655502,
                                    -631.301758441,
                                    0.15083915761],
                                   [u'p__Firmicutes',
                                    5.97534179688e-06,
                                    15,
                                    0.333322166804,
                                    -631.301771907,
                                    0.150839131805]]}],
                          u'tree': 'tree_thanks',
                          u'version': 3}

        with tempdir.TempDir() as d:
            cmd = "%s pipe --sequences %s --otu_table /dev/null --output_jplace %s"\
                  " --singlem_packages %s" % (
                      path_to_script,
                      os.path.join(path_to_data,'1_pipe','jplace_test.fna'),
                      os.path.join(d, "my_jplace"),
                      os.path.join(path_to_data,'4.12.22seqs.spkg'))
            extern.run(cmd)
            jplace_path = os.path.join(d, 'my_jplace_jplace_test_4.12.22seqs.jplace')
            j = json.load(open(jplace_path))
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
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_revcom_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	1	2.56	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
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
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
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
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s %s" % (
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
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s --no_assign_taxonomy" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_known_sequence_taxonomy(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	mytax; yeh',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>another
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile() as taxf:
                taxf.write("HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482\tmytax; yeh\n")
                taxf.write("another\tmytax; yeh; 2\n")
                taxf.flush()

                cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s "\
                      "--no_assign_taxonomy --known_sequence_taxonomy %s"% (
                          path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                          taxf.name)
                self.assertEqual(expected,
                                 extern.run(cmd).replace(
                                     os.path.basename(n.name).replace('.fa',''),
                                     '').split("\n"))

    def test_sample_name_strange_characters(self):
        expected = [self.headers,
                    ['4.12.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.8',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.8',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes'],
                    ['4.12.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.9',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.9',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes']]
        cmd = "%s --quiet pipe --sequences "\
        "%s/1_pipe/contigs.fasta.metabat-bins-_-t20_--superspecific.9.fa "\
        "%s/1_pipe/contigs.fasta.metabat-bins-_-t20_--superspecific.8.fa "\
        "--otu_table /dev/stdout --threads 4 --singlem_packages %s" %(
            path_to_script,
            path_to_data,
            path_to_data,
            self.two_packages)
        self.assertEqualOtuTable(expected, extern.run(cmd))

    def test_archive_otu_groopm_compatibility(self):
        """This tests for API stability, where the API is used by GroopM 2.0"""
        expected = [('contig_1', '4.11.22seqs', 'Root; d__Bacteria; p__Firmicutes')]

        inseqs = '''>contig_1 abcd
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --archive_otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))

            j = json.loads(extern.run(cmd))
            fields = j['fields']
            data = j['otus']
            self.assertEqual(expected,
                             [(name, row[fields.index('gene')], row[fields.index('taxonomy')]) for row in data for name in row[fields.index('read_names')]]
                            )

    def test__align_proteins_to_hmm(self):
        proteins = path_to_data+'/4.12.22seqs.spkg/4.12.22seqs/singlem_package_creatorq4droc.fasta'
        hmm = path_to_data+'/4.12.22seqs.spkg/4.12.22seqs/graftmgyqgXl_search.hmm'

        alignment = SearchPipe()._align_proteins_to_hmm(proteins, hmm)
        self.assertEqual(22, len(alignment))
        a = alignment[0]
        self.assertEqual('2512564006', a.name)
        self.assertEqual(
            '-------MAKKVAGTMKLQVAAGKANPSPPVGPALGQRGINIMEFCKAFNAKTaDLEP-----GAPCPTVITYYQDKSFSMEIKTPPASYFLKKAAKV-----K--------SGSKTPSRDTVG---------TVTTKQVREIAEAKMKDLNANDIEGAMKIILGSARSMGIEVK---------',
            a.seq)
        a2 = alignment[21]
        self.assertEqual('2519103189', a2.name)
        self.assertEqual(
            '-------VAKKVDSVVKLQIPAGKANPAPPVGPALGQAGINIMGFCKEFNAQT-QDQA-----GMIIPVEITVYEDRSFTFITKTPPAAVLLKKAAGI-----E--------TASGEPNRNKVA---------TLNRDKVKEIAELKMPDLNAADVEAAMRMVEGTARSMGIVIED--------',
            a2.seq)

    def test_protein_package_non60_length(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTT	1	1.41	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
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
        with tempfile.NamedTemporaryFile(suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.length10.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))


if __name__ == "__main__":
    unittest.main()
