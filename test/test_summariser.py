#!/usr/bin/env python

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

import unittest
import os.path
from string import split
import sys
from StringIO import StringIO
import tempfile
import tempdir
import extern

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.summariser import Summariser
from singlem.otu_table_collection import OtuTableCollection

class Tests(unittest.TestCase):
    headers = split('gene sample sequence num_hits coverage taxonomy')
    output_headers = split('type gene sample difference_in_bp sequence num_hits coverage taxonomy')

    def test_archive_to_otu_table_conversion(self):
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
        with tempdir.TempDir() as tmp:
            Summariser.summarise(krona_output=os.path.join(tmp, 'KronaOK.html'),
                                 table_collection=table_collection)
            self.assertTrue(os.path.exists(os.path.join(tmp,'KronaOK.html')))

    def test_biom_hello_world(self):
        insert_otu_table = [self.headers,
                            ['4.12.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                            ['4.12.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTtttCAAGCAGGTGTG','2','2.94','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        with tempdir.TempDir() as tmp:
            with tempfile.NamedTemporaryFile(suffix='.otu_table.csv') as n:
                n.write("\n".join(["\t".join(x) for x in insert_otu_table]+['']))
                n.flush()
                extern.run("%s summarise --biom_prefix '%s' --input_otu_tables '%s'" % (
                    path_to_script, os.path.join(tmp,"mybiom"), n.name))
                self.assertEqual(['mybiom.4.12.ribosomal_protein_L11_rplK.biom'], os.listdir(tmp))
                self.assertEqual(
                    '# Constructed from biom file\n#OTU ID\tinsert\ttaxonomy\nRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG\t1.0\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales\nRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTtttCAAGCAGGTGTG\t2.0\tRoot; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales',
                    extern.run("biom convert -i '%s' -o /dev/stdout --to-tsv --header-key taxonomy" % os.path.join(tmp,'mybiom.4.12.ribosomal_protein_L11_rplK.biom')))

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


if __name__ == "__main__":
    unittest.main()
