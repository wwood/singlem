#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft
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

import json
import unittest
import os.path
import sys
from io import StringIO
import tempfile
import extern
import re

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.summariser import Summariser
from singlem.otu_table_collection import OtuTableCollection
from singlem.archive_otu_table import ArchiveOtuTable

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    output_headers = str.split('type gene sample difference_in_bp sequence num_hits coverage taxonomy')
    maxDiff = None

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))

    def test_archive_to_otu_table_conversion(self):
        cmd = "{} summarise --input-gzip-archive-otu-table-list <(ls {}/small.otu_table.json.gz) --output-otu-table /dev/stdout".format(
            path_to_script,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy\n\
S1.5.ribosomal_protein_L11_rplK\tsmall\tCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG\t4\t9.76\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\n', obs)

    def test_gzip_archive_list_to_otu_table_conversion(self):

        archive = '{"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?"], "singlem_package_sha256s": ["2b2afe0114de20451fccfe74360756376dc83d001d890e84e322ab0833eca6ba", "7f406a73d8bb176994055cb966ff350e208986d12c8215722686c17c26e548c7", "735b44ae547c133163cb7d40f417292c35423864d00c95e7f1b32091b27d46c5", "8fc6dcce2766cc01defb3b5c689a1ed8ce9d59b725c67e58c2044dafaae908b3", "172df49937742b8411d41d217500d862567374401eaf393b25107b22ac630202", "4cb1bf226bf28d8198ed5c29e8a76df411d96a6c3ce1256af16887b9a184b0a6", "d473d3ae677e6e46202461ccdedb2aef23c0a10a3412422586b37e397ca37294", "431a2860bb890cd1c7193c565cbf0cc227850cba36fb17fe94df686e74ee9b11", "faa663527bb9aea63cef03859311f2e7f55fe98590a5ec85c5ba85815a6fd13e", "a0daf111380e6e499ad9c10c3ac413aa9016c7503dd459825100168524bff0d1", "aba631d4735aeb9d2dfbbbfee1c0739bf9e99ad6532a3be04ff627f3e6efdae2", "bba10c1feb0c26bdf46aa3d1dcb992744a699cde5cf02bb2728f8397378b342f", "4d91dd794b25fd256508f0814f6a2d31e20dc85e0aa9ea405031398565276768", "9b23c524a6210af0706eea7252c2d378888029f141b9305c3e88cbac3fd83f88", "50a209417b455a48bc67702d6a3809a172c57f00785d8c705a3322e6e2a71f72"], "version": 1, "alignment_hmm_sha256s": ["dd9b7e283598360b89ec91ff3f5c509361a6108a2eadc44bfb29646b1510f6b7", "b1bb943c3449a78f937db960bfdf6b2bed641388d33fce3cb2d5f69e79946ea6", "de92c90f2c83e380ae3953972fb63fcb8ce868dab87a305f9f1811b84ffb3d39", "453ed4a62608a4aec36117a2dd1a276709ff6b130ecb8d7b1612926bfab25527", "20cc450cf4157ecf1772e0325d4d8ed400b597d888a5cb5044ca69098f935656", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "7cbba7ba0ed58d21c7519ba3fcef0abe43378e5c38c985b0d5e0e5219f141d92", "4a3bbe5ac594ef3c7c820e74544828e19eca68bf860d64f928729eb4530fce4e", "06a4bed0a765971b891ca4a4bf5680aeef4a4a249ce0c028798c0e912f0ccfb4", "2678fe218ca860a2d88bdbf76935d8c78a00ab6603a041a432505d754ef08250", "b54ff98aa03ab31af39c737a569b23ee4ed9296c8ea088562bfb3db87c38fe4a", "4ae31f14067bf183f38dca20f2aefb580e5ff25848881dd988908b70b67761bb", "d7bb3d544133f38110a329712b3ace7e7d7c989dafa3815d2d5a292b4c575f50", "7639bb919ef54f7baff3ed3a8c924efca97ed375cf4120a6e05d98fd6ef52cbb", "6923b889888ea34fabf463b2c8ad5fe23c94828f1a2631a07601f246f5e87150"], "otus": [["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false], ["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false]]}'
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.add_archive_otu_table(StringIO(archive))
        Summariser().write_otu_table(
            table_collection = table_collection,
            output_table_io = output,
            output_extras=False)
        self.assertEqual(exp, output.getvalue())

    def test_archive_to_otu_table_output_extras(self):
        archive = '{"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?"], "singlem_package_sha256s": ["2b2afe0114de20451fccfe74360756376dc83d001d890e84e322ab0833eca6ba", "7f406a73d8bb176994055cb966ff350e208986d12c8215722686c17c26e548c7", "735b44ae547c133163cb7d40f417292c35423864d00c95e7f1b32091b27d46c5", "8fc6dcce2766cc01defb3b5c689a1ed8ce9d59b725c67e58c2044dafaae908b3", "172df49937742b8411d41d217500d862567374401eaf393b25107b22ac630202", "4cb1bf226bf28d8198ed5c29e8a76df411d96a6c3ce1256af16887b9a184b0a6", "d473d3ae677e6e46202461ccdedb2aef23c0a10a3412422586b37e397ca37294", "431a2860bb890cd1c7193c565cbf0cc227850cba36fb17fe94df686e74ee9b11", "faa663527bb9aea63cef03859311f2e7f55fe98590a5ec85c5ba85815a6fd13e", "a0daf111380e6e499ad9c10c3ac413aa9016c7503dd459825100168524bff0d1", "aba631d4735aeb9d2dfbbbfee1c0739bf9e99ad6532a3be04ff627f3e6efdae2", "bba10c1feb0c26bdf46aa3d1dcb992744a699cde5cf02bb2728f8397378b342f", "4d91dd794b25fd256508f0814f6a2d31e20dc85e0aa9ea405031398565276768", "9b23c524a6210af0706eea7252c2d378888029f141b9305c3e88cbac3fd83f88", "50a209417b455a48bc67702d6a3809a172c57f00785d8c705a3322e6e2a71f72"], "version": 1, "alignment_hmm_sha256s": ["dd9b7e283598360b89ec91ff3f5c509361a6108a2eadc44bfb29646b1510f6b7", "b1bb943c3449a78f937db960bfdf6b2bed641388d33fce3cb2d5f69e79946ea6", "de92c90f2c83e380ae3953972fb63fcb8ce868dab87a305f9f1811b84ffb3d39", "453ed4a62608a4aec36117a2dd1a276709ff6b130ecb8d7b1612926bfab25527", "20cc450cf4157ecf1772e0325d4d8ed400b597d888a5cb5044ca69098f935656", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "7cbba7ba0ed58d21c7519ba3fcef0abe43378e5c38c985b0d5e0e5219f141d92", "4a3bbe5ac594ef3c7c820e74544828e19eca68bf860d64f928729eb4530fce4e", "06a4bed0a765971b891ca4a4bf5680aeef4a4a249ce0c028798c0e912f0ccfb4", "2678fe218ca860a2d88bdbf76935d8c78a00ab6603a041a432505d754ef08250", "b54ff98aa03ab31af39c737a569b23ee4ed9296c8ea088562bfb3db87c38fe4a", "4ae31f14067bf183f38dca20f2aefb580e5ff25848881dd988908b70b67761bb", "d7bb3d544133f38110a329712b3ace7e7d7c989dafa3815d2d5a292b4c575f50", "7639bb919ef54f7baff3ed3a8c924efca97ed375cf4120a6e05d98fd6ef52cbb", "6923b889888ea34fabf463b2c8ad5fe23c94828f1a2631a07601f246f5e87150"], "otus": [["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false], ["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false]]}'
        exp = 'gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy\tread_names\tnucleotides_aligned\ttaxonomy_by_known?\n4.11.ribosomal_protein_L10\tminimal\tTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA\t2\t4.88\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus\tHWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331\t60 60\tFalse\n4.12.ribosomal_protein_L11_rplK\tminimal\tCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG\t4\t9.76\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\tHWI-ST1243:156:D1K83ACXX:7:1109:18214:9910 HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124 HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928 HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842\t60 60 60 60\tFalse\n'
        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.add_archive_otu_table(StringIO(archive))
        Summariser().write_otu_table(
            table_collection = table_collection,
            output_table_io = output,
            output_extras=True)
        self.assertEqual(exp, output.getvalue())

    def test_two_archives_to_otu_table(self):
        archive = '{"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?"], "singlem_package_sha256s": ["2b2afe0114de20451fccfe74360756376dc83d001d890e84e322ab0833eca6ba", "7f406a73d8bb176994055cb966ff350e208986d12c8215722686c17c26e548c7", "735b44ae547c133163cb7d40f417292c35423864d00c95e7f1b32091b27d46c5", "8fc6dcce2766cc01defb3b5c689a1ed8ce9d59b725c67e58c2044dafaae908b3", "172df49937742b8411d41d217500d862567374401eaf393b25107b22ac630202", "4cb1bf226bf28d8198ed5c29e8a76df411d96a6c3ce1256af16887b9a184b0a6", "d473d3ae677e6e46202461ccdedb2aef23c0a10a3412422586b37e397ca37294", "431a2860bb890cd1c7193c565cbf0cc227850cba36fb17fe94df686e74ee9b11", "faa663527bb9aea63cef03859311f2e7f55fe98590a5ec85c5ba85815a6fd13e", "a0daf111380e6e499ad9c10c3ac413aa9016c7503dd459825100168524bff0d1", "aba631d4735aeb9d2dfbbbfee1c0739bf9e99ad6532a3be04ff627f3e6efdae2", "bba10c1feb0c26bdf46aa3d1dcb992744a699cde5cf02bb2728f8397378b342f", "4d91dd794b25fd256508f0814f6a2d31e20dc85e0aa9ea405031398565276768", "9b23c524a6210af0706eea7252c2d378888029f141b9305c3e88cbac3fd83f88", "50a209417b455a48bc67702d6a3809a172c57f00785d8c705a3322e6e2a71f72"], "version": 1, "alignment_hmm_sha256s": ["dd9b7e283598360b89ec91ff3f5c509361a6108a2eadc44bfb29646b1510f6b7", "b1bb943c3449a78f937db960bfdf6b2bed641388d33fce3cb2d5f69e79946ea6", "de92c90f2c83e380ae3953972fb63fcb8ce868dab87a305f9f1811b84ffb3d39", "453ed4a62608a4aec36117a2dd1a276709ff6b130ecb8d7b1612926bfab25527", "20cc450cf4157ecf1772e0325d4d8ed400b597d888a5cb5044ca69098f935656", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "7cbba7ba0ed58d21c7519ba3fcef0abe43378e5c38c985b0d5e0e5219f141d92", "4a3bbe5ac594ef3c7c820e74544828e19eca68bf860d64f928729eb4530fce4e", "06a4bed0a765971b891ca4a4bf5680aeef4a4a249ce0c028798c0e912f0ccfb4", "2678fe218ca860a2d88bdbf76935d8c78a00ab6603a041a432505d754ef08250", "b54ff98aa03ab31af39c737a569b23ee4ed9296c8ea088562bfb3db87c38fe4a", "4ae31f14067bf183f38dca20f2aefb580e5ff25848881dd988908b70b67761bb", "d7bb3d544133f38110a329712b3ace7e7d7c989dafa3815d2d5a292b4c575f50", "7639bb919ef54f7baff3ed3a8c924efca97ed375cf4120a6e05d98fd6ef52cbb", "6923b889888ea34fabf463b2c8ad5fe23c94828f1a2631a07601f246f5e87150"], "otus": [["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false], ["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false]]}'
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','maximal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','maximal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']
            ]
        exp = "\n".join(["\t".join(x) for x in e]+[''])

        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.add_archive_otu_table(StringIO(archive))
        table_collection.add_archive_otu_table(StringIO(archive.replace('minimal', 'maximal')))
        Summariser().write_otu_table(
            table_collection = table_collection,
            output_table_io = output,
            output_extras=False)
        self.assertEqual(exp, output.getvalue())

    def test_krona(self):
        archive = '{"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?"], "singlem_package_sha256s": ["2b2afe0114de20451fccfe74360756376dc83d001d890e84e322ab0833eca6ba", "7f406a73d8bb176994055cb966ff350e208986d12c8215722686c17c26e548c7", "735b44ae547c133163cb7d40f417292c35423864d00c95e7f1b32091b27d46c5", "8fc6dcce2766cc01defb3b5c689a1ed8ce9d59b725c67e58c2044dafaae908b3", "172df49937742b8411d41d217500d862567374401eaf393b25107b22ac630202", "4cb1bf226bf28d8198ed5c29e8a76df411d96a6c3ce1256af16887b9a184b0a6", "d473d3ae677e6e46202461ccdedb2aef23c0a10a3412422586b37e397ca37294", "431a2860bb890cd1c7193c565cbf0cc227850cba36fb17fe94df686e74ee9b11", "faa663527bb9aea63cef03859311f2e7f55fe98590a5ec85c5ba85815a6fd13e", "a0daf111380e6e499ad9c10c3ac413aa9016c7503dd459825100168524bff0d1", "aba631d4735aeb9d2dfbbbfee1c0739bf9e99ad6532a3be04ff627f3e6efdae2", "bba10c1feb0c26bdf46aa3d1dcb992744a699cde5cf02bb2728f8397378b342f", "4d91dd794b25fd256508f0814f6a2d31e20dc85e0aa9ea405031398565276768", "9b23c524a6210af0706eea7252c2d378888029f141b9305c3e88cbac3fd83f88", "50a209417b455a48bc67702d6a3809a172c57f00785d8c705a3322e6e2a71f72"], "version": 1, "alignment_hmm_sha256s": ["dd9b7e283598360b89ec91ff3f5c509361a6108a2eadc44bfb29646b1510f6b7", "b1bb943c3449a78f937db960bfdf6b2bed641388d33fce3cb2d5f69e79946ea6", "de92c90f2c83e380ae3953972fb63fcb8ce868dab87a305f9f1811b84ffb3d39", "453ed4a62608a4aec36117a2dd1a276709ff6b130ecb8d7b1612926bfab25527", "20cc450cf4157ecf1772e0325d4d8ed400b597d888a5cb5044ca69098f935656", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "7cbba7ba0ed58d21c7519ba3fcef0abe43378e5c38c985b0d5e0e5219f141d92", "4a3bbe5ac594ef3c7c820e74544828e19eca68bf860d64f928729eb4530fce4e", "06a4bed0a765971b891ca4a4bf5680aeef4a4a249ce0c028798c0e912f0ccfb4", "2678fe218ca860a2d88bdbf76935d8c78a00ab6603a041a432505d754ef08250", "b54ff98aa03ab31af39c737a569b23ee4ed9296c8ea088562bfb3db87c38fe4a", "4ae31f14067bf183f38dca20f2aefb580e5ff25848881dd988908b70b67761bb", "d7bb3d544133f38110a329712b3ace7e7d7c989dafa3815d2d5a292b4c575f50", "7639bb919ef54f7baff3ed3a8c924efca97ed375cf4120a6e05d98fd6ef52cbb", "6923b889888ea34fabf463b2c8ad5fe23c94828f1a2631a07601f246f5e87150"], "otus": [["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false], ["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false]]}'
        table_collection = OtuTableCollection()
        table_collection.add_archive_otu_table(StringIO(archive))
        with tempfile.TemporaryDirectory() as tmp:
            Summariser.write_otu_table_krona(krona_output=os.path.join(tmp, 'KronaOK.html'),
                                 table_collection=table_collection)
            self.assertTrue(os.path.exists(os.path.join(tmp,'KronaOK.html')))

    def test_wide_format(self):
        e = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.11.ribosomal_protein_L10','minimal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','maximal','TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA','2','4.88','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        exp = "\n".join(["\t".join(x) for x in e]+[''])
        output = StringIO()
        table_collection = OtuTableCollection()
        table_collection.add_otu_table(StringIO(exp))
        Summariser().write_wide_format_otu_table(
            table_collection = table_collection,
            output_table_io = output)
        self.assertEqual('marker\tsequence\tminimal\tmaximal\ttaxonomy\n4.11.ribosomal_protein_L10\tTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA\t2\t2\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus\n4.12.ribosomal_protein_L11_rplK\tCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG\t4\t0\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\n',
                         output.getvalue())

    def test_phylogeny_aware_beta_diversity(self):
        '''Need to make sure we can do phylogeny aware stuff with the output of summarise as per the README'''
        e1 = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2512564006'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','ACTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2585427686'],
            ['4.12.ribosomal_protein_L11_rplK','minimal','TCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2518645625']]
        e2 = [['gene','sample','sequence','num_hits','coverage','taxonomy'],
            ['4.12.ribosomal_protein_L11_rplK','minimal2','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2512564006'],
            ['4.12.ribosomal_protein_L11_rplK','minimal2','ACTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2585427686']]

        with tempfile.NamedTemporaryFile(mode='w') as otus1:
            otus1.write("\n".join(["\t".join(x) for x in e1]+['']))
            otus1.flush()
            with tempfile.NamedTemporaryFile(mode='w') as otus2:
                otus2.write("\n".join(["\t".join(x) for x in e2]+['']))
                otus2.flush()

                with tempfile.TemporaryDirectory() as d:
                    #singlem summarise --input-otu-table /tmp/tmpGTvm6f /tmp/tmpTvoNvC --unifrac out_unifracer
                    cmd = "{} summarise --input-otu-tables {} {} --unifrac-by-otu {}".format(
                        path_to_script,
                        otus1.name,
                        otus2.name,
                        os.path.join(d, 'unifrac_out'))
                    extern.run(cmd)
                    #convertToEBD.py out_unifracer.4.12.ribosomal_protein_L11_rplK.unifrac otu_table.ebd
                    cmd = "convertToEBD.py {} {}".format(
                        os.path.join(d,"unifrac_out.4.12.ribosomal_protein_L11_rplK.unifrac"),
                        os.path.join(d, "otu_table1.ebd"))
                    extern.run(cmd)
                    #ExpressBetaDiversity -s otu_table.ebd -c Bray-Curtis
                    cmd = "ExpressBetaDiversity -s {} -c Bray-Curtis -p {}".format(
                        os.path.join(d, "otu_table1.ebd"),
                        os.path.join(d, "phylogeny_free"))
                    extern.run(cmd)
                    with open(os.path.join(d, "phylogeny_free.diss")) as f:
                        r = f.read()
                        self.assertTrue(r=="""2
minimal2
minimal	0.2
"""
or r=="""2
minimal
minimal2	0.2
""")

                with tempfile.TemporaryDirectory() as d:
                    ## Then for the phylogeny version:
                    #singlem summarise --taxonomy_as_identifier --unifrac
                    cmd = "{} summarise --input-otu-tables {} {} --unifrac-by-taxonomy {}".format(
                        path_to_script,
                        otus1.name,
                        otus2.name,
                        os.path.join(d, 'unifrac_taxonomy_out'))
                    extern.run(cmd)
                    cmd = "convertToEBD.py {} {}".format(
                        os.path.join(d,"unifrac_taxonomy_out.4.12.ribosomal_protein_L11_rplK.unifrac"),
                        os.path.join(d, "otu_table2.ebd"))
                    extern.run(cmd)
                    cmd = "ExpressBetaDiversity -s {} -c Bray-Curtis -t {} -p {}".format(
                        os.path.join(d, 'otu_table2.ebd'),
                        os.path.join(
                            path_to_data,
                            "4.12.22seqs.spkg/4.12.22seqs/4.12.22seqs.gpkg.refpkg/treeitF_La.tre"),
                        os.path.join(d, 'phylogeny_full'))
                    extern.run(cmd)
                    with open(os.path.join(d, "phylogeny_full.diss")) as f:
                        r = f.read()
                        self.assertTrue(r=='2\nminimal2\nminimal\t0.109937\n' or r=='2\nminimal\nminimal2\t0.109937\n')

    @unittest.skip("Default packages no longer in git repo, and may not contain trees anyway, unsure")
    def test_get_tree_default(self):
        cmd = "{} get_tree".format(path_to_script)
        observed = extern.run(cmd)
        splits = observed.split('\n')
        self.assertEqual('marker\ttree_file', splits[0])
        self.assertEqual('.tre',splits[1][-4:])
        self.assertGreater(len(splits), 10)
        for line in splits[1:-1]:
            self.assertTrue(os.path.exists(line.split('\t')[1]))

    def test_get_tree_specific(self):
        cmd = "{} get_tree --singlem-package {}".format(
            path_to_script,
            os.path.join(path_to_data,'4.12.22seqs.spkg'))
        observed = extern.run(cmd)
        splits = observed.split('\n')
        self.assertEqual('marker\ttree_file', splits[0])
        self.assertEqual(3, len(splits))
        self.assertEqual('',splits[2])
        self.assertEqual(
            '4.12.22seqs\t{}'.format(
                os.path.abspath(os.path.join(
                    path_to_data,
                    '4.12.22seqs.spkg/4.12.22seqs/4.12.22seqs.gpkg.refpkg/treeitF_La.tre'))),
            splits[1])

    def test_exclude_off_target_domains_keep(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs	small	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes',
            '']
        cmd = "{} summarise --singlem-packages {} --exclude-off-target-hits --input-archive-otu-table {} --output-otu-table /dev/stdout".format(
            path_to_script,
            os.path.join(path_to_data,'4.11.22seqs.v3.gpkg.spkg'),
            os.path.join(path_to_data,'small.otu_table.4.11.22seqs.json'))
        observed = extern.run(cmd)
        self.assertEqualOtuTable(list([line.split("\t") for line in expected]), observed)

    def test_exclude_off_target_domains_exclude(self):
        expected = [
            "\t".join(self.headers),
            '']
        cmd = "{} summarise --singlem-packages {} --exclude-off-target-hits --input-archive-otu-table {} --output-otu-table /dev/stdout".format(
            path_to_script,
            os.path.join(path_to_data,'4.11.22seqs.v3_archaea_targetted.gpkg.spkg'),
            os.path.join(path_to_data,'small.otu_table.4.11.22seqs.json'))
        observed = extern.run(cmd)
        self.assertEqualOtuTable(list([line.split("\t") for line in expected]), observed)

    def test_translate(self):
        expected = [
            "\t".join(self.headers),
            '4.11.ribosomal_protein_L10	small	LRSQLREAGVEYKVYKNTMV	4	9.76	Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus',
            '4.11.ribosomal_protein_L10	small	KRSQLREAGVEYKVYKNTMV	5	10.4	Root; d__Bacteria',
            '4.12.ribosomal_protein_L11_rplK	small	PAGKANPAPPVGPALGQAGV	4	9.76	Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales',
            '']
        cmd = "{} summarise --input-otu-table {} --output-translated-otu-table /dev/stdout".format(
            path_to_script,
            os.path.join(path_to_data,'small_will_collapse.otu_table.csv'))
        observed = extern.run(cmd)
        self.assertEqualOtuTable(list([line.split("\t") for line in expected]), observed)

    def test_collapse_paired_with_unpaired(self):
        cmd = '{} summarise --input-archive-otu-table {}/summarise/ERR1883421_paired.annotated.singlem.renew.renew_v4_renew.json {}/summarise/ERR1883421_unpaired.annotated.singlem.renew.renew_v4_renew.json --collapse-paired-with-unpaired /dev/stdout'.format(
            path_to_script,
            path_to_data,
            path_to_data)
        res = extern.run(cmd)
        ar = ArchiveOtuTable.read(StringIO(res))
        self.assertEqual(len(ar.data), 565)
        
    def test_dump_raw_sequences(self):
        cmd = '{} summarise --input-archive-otu-table {}/small.otu_table.4.11.22seqs.json --unaligned-sequences-dump-file /dev/stdout'.format(
            path_to_script,
            path_to_data)
        res = extern.run(cmd)
        expected = """
            >HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331~0
            ACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAATCCACGGTAGTCA
            >HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482~0
            ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
            """
        expected = '\n'.join([line.strip() for line in expected.split('\n') if line.strip()])+"\n"
        self.assertEqual(expected, res)

    def test_collapse_to_sample_name(self):
        with tempfile.NamedTemporaryFile() as tf:
            extern.run(f'bin/singlem summarise --collapse-to-sample-name a1 --input-archive-otu-tables {path_to_data}/summarise/original1_for_merge.json  {path_to_data}/summarise/original2_for_merge.json --output-archive-otu-table {tf.name}')

            with open(tf.name) as f:
                observed = json.load(f)
                with open(f'{path_to_data}/summarise/collapsed_to_sample_name.json') as f:
                    expected = json.load(f)
                    self.assertEqual(observed, expected)

    def test_collapse_to_sample_name_no_assign_taxonomy(self):
        with tempfile.NamedTemporaryFile() as tf:
            extern.run(f'bin/singlem summarise --collapse-to-sample-name testsample --input-archive-otu-tables {path_to_data}/small.otu_table.no_assign_taxonomy.json  {path_to_data}/small.otu_table.no_assign_taxonomy.json --output-archive-otu-table {tf.name}')

            with open(tf.name) as f:
                observed = json.load(f)
                with open(f'{path_to_data}/summarise/collapsed_to_sample_name.no_assign_taxonomy.json') as f:
                    expected = json.load(f)
                    self.assertEqual(observed, expected)

    def test_species_by_site_relative(self):
        stdout = extern.run(f'bin/singlem summarise --input-taxonomic-profile {path_to_data}/summarise/marine0.head5.profile \
            --output-species-by-site-relative-abundance /dev/stdout \
            --output-species-by-site-level phylum')
        self.assertEqual(stdout, """taxonomy\tmarine0.1
unassigned\t50.77
Root; d__Archaea; p__Thermoproteota\t7.81
Root; d__Bacteria; p__Desulfobacterota\t11.16
Root; d__Bacteria; p__Proteobacteria\t30.26
""")

    def test_species_by_site_2_samples(self):
        stdout = extern.run(f'bin/singlem summarise --input-taxonomic-profile {path_to_data}/summarise/marine0.head5.profile \
             {path_to_data}/summarise/land.profile \
            --output-species-by-site-relative-abundance /dev/stdout \
            --output-species-by-site-level phylum')

        expected = re.compile(r'  +').sub('\t', """taxonomy        marine0.1       land0.1
unassigned      50.77   0.0
Root; d__Archaea; p__Thermoproteota     7.81    0.0
Root; d__Bacteria; p__Desulfobacterota  11.16   0.0
Root; d__Bacteria; p__Proteobacteria    30.26   14.58
Root; d__Bacteria; p__Firmicutes        0.0     75.61
Root; d__Bacteria; p__Firmicutes_A      0.0     6.98
Root; d__Bacteria; p__Actinobacteriota  0.0     2.83
""")
        self.assertEqual(expected, stdout)

    def test_concatenate_profiles(self):
        stdout = extern.run(f'bin/singlem summarise --input-taxonomic-profile {path_to_data}/summarise/profile1.tsv \
             {path_to_data}/summarise/profile2.tsv \
            --output-taxonomic-profile /dev/stdout')

        expected = re.compile(r'  +').sub('\t', """sample        coverage       taxonomy
A2A_APEX_ALPHA_Plot_21_A30_10.1  18.69   Root; d__Archaea
A2A_APEX_ALPHA_Plot_21_A30_10.1  89.04   Root; d__Bacteria
A2A_APEX_ALPHA_Plot_21_A30_10.1  3.57    Root; d__Archaea; p__Halobacteriota
sample2  14.0      Root; d__Archaea
sample2  15.0     Root; d__Bacteria
""")
        self.assertEqual(expected, stdout)

    def test_concatenate_profiles_duplicate_samples(self):
        with self.assertRaises(Exception) as cm:
            extern.run(f'bin/singlem summarise --input-taxonomic-profile {path_to_data}/summarise/profile1.tsv \
             {path_to_data}/summarise/profile1.tsv \
            --output-taxonomic-profile /dev/stdout')

    def test_fill_condensed(self):
        stdout = extern.run(f'bin/singlem summarise --input-taxonomic-profile <(head -5 {path_to_data}/read_fraction/marine0.profile) <(head -5 {path_to_data}/read_fraction/marine0.profile |sed s/marine0.1/land/) '+
            '--output-filled-taxonomic-profile /dev/stdout')

        expected = re.compile(r'  +').sub('\t', """sample  filled_coverage        taxonomy
marine0.1       7.17    Root
marine0.1       4.20    Root; d__Archaea
marine0.1       2.97    Root; d__Bacteria
marine0.1       0.56    Root; d__Archaea; p__Thermoproteota
marine0.1       0.80     Root; d__Bacteria; p__Desulfobacterota
marine0.1       2.17    Root; d__Bacteria; p__Proteobacteria
land       7.17    Root
land       4.20    Root; d__Archaea
land       2.97    Root; d__Bacteria
land       0.56    Root; d__Archaea; p__Thermoproteota
land       0.80     Root; d__Bacteria; p__Desulfobacterota
land       2.17    Root; d__Bacteria; p__Proteobacteria
""")
        self.assertEqual(expected, stdout)
        

if __name__ == "__main__":
    unittest.main()
