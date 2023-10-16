
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

import unittest
import os.path
import sys
from io import StringIO
import tempfile
import extern
from bird_tool_utils import in_tempdir

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
DEFAULT_WINDOW_SIZE = 60

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.appraiser import Appraiser
from singlem.otu_table_collection import OtuTableCollection
from singlem.otu_table import OtuTable
from singlem.metapackage import Metapackage

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None

    def setUp(self):
        OtuTable._clear_cache()

    def _sort_appraisal_results(self, appraisal_results):
        sorted_appraisal_results = list(sorted(
            appraisal_results,
            key=lambda x: (x.metagenome_sample_name, x.num_binned, x.num_assembled, x.num_not_found)))
        return sorted_appraisal_results

    def test_hello_world(self):
        metagenome_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(1, len(a.binned_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.binned_otus[0].sequence)
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                         a.not_found_otus[0].sequence)

    def test_hello_word_cmdline(self):
        with in_tempdir():
            cmd = (
                "{} appraise "
                "--metagenome-otu-tables {}/appraise_example3/reads.otu_table.tsv "
                "--genome-otu-tables {}/appraise_example3/bins.otu_table.tsv "
                "--metapackage {}/four_package.smpkg "
                "--output-binned-otu-table binned.otu_table.tsv "
                "--output-unaccounted-for-otu-table unbinned.otu_table.tsv "
            ).format(
                path_to_script,
                path_to_data,
                path_to_data,
                path_to_data
            )
            extern.run(cmd)

            expected_binned = [
                self.headers,
                ["4.12.ribosomal_protein_L11_rplK", "minimal", "GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC", "7", "17.07", "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales"],
                ""
            ]
            expected_binned = "\n".join(["\t".join(x) for x in expected_binned])

            expected_unbinned = [
                self.headers,
                ["4.11.ribosomal_protein_L10", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", "4", "9.76", "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus"],
                ""
            ]
            expected_unbinned = "\n".join(["\t".join(x) for x in expected_unbinned])

            with open('binned.otu_table.tsv') as f:
                observed_binned = "".join(f.readlines())
            with open('unbinned.otu_table.tsv') as f:
                observed_unbinned = "".join(f.readlines())

            self.assertEqual(expected_binned, observed_binned)
            self.assertEqual(expected_unbinned, observed_unbinned)

    def test_multiple_samples(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)
        a = res[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)

    def test_archive_input(self):
        metagenomes = '{"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?"], "singlem_package_sha256s": ["2b2afe0114de20451fccfe74360756376dc83d001d890e84e322ab0833eca6ba", "7f406a73d8bb176994055cb966ff350e208986d12c8215722686c17c26e548c7", "735b44ae547c133163cb7d40f417292c35423864d00c95e7f1b32091b27d46c5", "8fc6dcce2766cc01defb3b5c689a1ed8ce9d59b725c67e58c2044dafaae908b3", "172df49937742b8411d41d217500d862567374401eaf393b25107b22ac630202", "4cb1bf226bf28d8198ed5c29e8a76df411d96a6c3ce1256af16887b9a184b0a6", "d473d3ae677e6e46202461ccdedb2aef23c0a10a3412422586b37e397ca37294", "431a2860bb890cd1c7193c565cbf0cc227850cba36fb17fe94df686e74ee9b11", "faa663527bb9aea63cef03859311f2e7f55fe98590a5ec85c5ba85815a6fd13e", "a0daf111380e6e499ad9c10c3ac413aa9016c7503dd459825100168524bff0d1", "aba631d4735aeb9d2dfbbbfee1c0739bf9e99ad6532a3be04ff627f3e6efdae2", "bba10c1feb0c26bdf46aa3d1dcb992744a699cde5cf02bb2728f8397378b342f", "4d91dd794b25fd256508f0814f6a2d31e20dc85e0aa9ea405031398565276768", "9b23c524a6210af0706eea7252c2d378888029f141b9305c3e88cbac3fd83f88", "50a209417b455a48bc67702d6a3809a172c57f00785d8c705a3322e6e2a71f72"], "version": 1, "alignment_hmm_sha256s": ["dd9b7e283598360b89ec91ff3f5c509361a6108a2eadc44bfb29646b1510f6b7", "b1bb943c3449a78f937db960bfdf6b2bed641388d33fce3cb2d5f69e79946ea6", "de92c90f2c83e380ae3953972fb63fcb8ce868dab87a305f9f1811b84ffb3d39", "453ed4a62608a4aec36117a2dd1a276709ff6b130ecb8d7b1612926bfab25527", "20cc450cf4157ecf1772e0325d4d8ed400b597d888a5cb5044ca69098f935656", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "7cbba7ba0ed58d21c7519ba3fcef0abe43378e5c38c985b0d5e0e5219f141d92", "4a3bbe5ac594ef3c7c820e74544828e19eca68bf860d64f928729eb4530fce4e", "06a4bed0a765971b891ca4a4bf5680aeef4a4a249ce0c028798c0e912f0ccfb4", "2678fe218ca860a2d88bdbf76935d8c78a00ab6603a041a432505d754ef08250", "b54ff98aa03ab31af39c737a569b23ee4ed9296c8ea088562bfb3db87c38fe4a", "4ae31f14067bf183f38dca20f2aefb580e5ff25848881dd988908b70b67761bb", "d7bb3d544133f38110a329712b3ace7e7d7c989dafa3815d2d5a292b4c575f50", "7639bb919ef54f7baff3ed3a8c924efca97ed375cf4120a6e05d98fd6ef52cbb", "6923b889888ea34fabf463b2c8ad5fe23c94828f1a2631a07601f246f5e87150"], "otus": [["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false], ["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false]]}'

        genomes = '{"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?"], "singlem_package_sha256s": ["2b2afe0114de20451fccfe74360756376dc83d001d890e84e322ab0833eca6ba", "7f406a73d8bb176994055cb966ff350e208986d12c8215722686c17c26e548c7", "735b44ae547c133163cb7d40f417292c35423864d00c95e7f1b32091b27d46c5", "8fc6dcce2766cc01defb3b5c689a1ed8ce9d59b725c67e58c2044dafaae908b3", "172df49937742b8411d41d217500d862567374401eaf393b25107b22ac630202", "4cb1bf226bf28d8198ed5c29e8a76df411d96a6c3ce1256af16887b9a184b0a6", "d473d3ae677e6e46202461ccdedb2aef23c0a10a3412422586b37e397ca37294", "431a2860bb890cd1c7193c565cbf0cc227850cba36fb17fe94df686e74ee9b11", "faa663527bb9aea63cef03859311f2e7f55fe98590a5ec85c5ba85815a6fd13e", "a0daf111380e6e499ad9c10c3ac413aa9016c7503dd459825100168524bff0d1", "aba631d4735aeb9d2dfbbbfee1c0739bf9e99ad6532a3be04ff627f3e6efdae2", "bba10c1feb0c26bdf46aa3d1dcb992744a699cde5cf02bb2728f8397378b342f", "4d91dd794b25fd256508f0814f6a2d31e20dc85e0aa9ea405031398565276768", "9b23c524a6210af0706eea7252c2d378888029f141b9305c3e88cbac3fd83f88", "50a209417b455a48bc67702d6a3809a172c57f00785d8c705a3322e6e2a71f72"], "version": 1, "alignment_hmm_sha256s": ["dd9b7e283598360b89ec91ff3f5c509361a6108a2eadc44bfb29646b1510f6b7", "b1bb943c3449a78f937db960bfdf6b2bed641388d33fce3cb2d5f69e79946ea6", "de92c90f2c83e380ae3953972fb63fcb8ce868dab87a305f9f1811b84ffb3d39", "453ed4a62608a4aec36117a2dd1a276709ff6b130ecb8d7b1612926bfab25527", "20cc450cf4157ecf1772e0325d4d8ed400b597d888a5cb5044ca69098f935656", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "7cbba7ba0ed58d21c7519ba3fcef0abe43378e5c38c985b0d5e0e5219f141d92", "4a3bbe5ac594ef3c7c820e74544828e19eca68bf860d64f928729eb4530fce4e", "06a4bed0a765971b891ca4a4bf5680aeef4a4a249ce0c028798c0e912f0ccfb4", "2678fe218ca860a2d88bdbf76935d8c78a00ab6603a041a432505d754ef08250", "b54ff98aa03ab31af39c737a569b23ee4ed9296c8ea088562bfb3db87c38fe4a", "4ae31f14067bf183f38dca20f2aefb580e5ff25848881dd988908b70b67761bb", "d7bb3d544133f38110a329712b3ace7e7d7c989dafa3815d2d5a292b4c575f50", "7639bb919ef54f7baff3ed3a8c924efca97ed375cf4120a6e05d98fd6ef52cbb", "6923b889888ea34fabf463b2c8ad5fe23c94828f1a2631a07601f246f5e87150"], "otus": [["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false]]}'

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_archive_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_archive_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[0]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(4.88/4)}, a.num_not_found)

    def test_clusterer_all_cluster_all_good(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATT','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.5,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(26.83/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)


    def test_clusterer_all_cluster_some_good(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','minimal',     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.7,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)


    def test_clusterer_all_cluster_two_samples(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','maximal',     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.7,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        res = sorted(app.appraisal_results)
        self.assertEqual(2, len(app.appraisal_results))
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)
        self.assertEqual(1, len(a.binned_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.binned_otus[0].sequence)
        self.assertEqual(0, len(a.not_found_otus))
        a = res[0]
        self.assertEqual('maximal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                         a.not_found_otus[0].sequence)
        self.assertEqual(0, len(a.binned_otus))


    def test_clusterer_all_cluster_two_samples_some_cluster(self):
        # non-As and genome cluster together but are not exactly the same
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.12.ribosomal_protein_L11_rplK','minimal','AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','12','24.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','maximal',     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA','4','8.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
                    ['4.11.ribosomal_protein_L10','maximal',     'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATG','5','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATA','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
                    ['4.11.ribosomal_protein_L10','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATT','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.7,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(24.07/4)}, a.num_not_found)
        self.assertEqual(1, len(a.binned_otus))
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.binned_otus[0].sequence)
        self.assertEqual('minimal',
                         a.binned_otus[0].sample_name)
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                         a.not_found_otus[0].sequence)
        self.assertEqual('minimal',
                         a.not_found_otus[0].sample_name)

        a = res[0]
        self.assertEqual('maximal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(8.76/4)}, a.num_not_found)
        self.assertEqual(1, len(a.binned_otus))
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATG',
                         a.binned_otus[0].sequence)
        self.assertEqual('maximal',
                         a.binned_otus[0].sample_name)
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                         a.not_found_otus[0].sequence)
        self.assertEqual('maximal',
                         a.not_found_otus[0].sample_name)

    def test_print_appraisal(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)
        a = res[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)

        to_print = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print)
        self.assertEqual("sample\tdomain\tnum_binned\tnum_not_found\tpercent_binned\nminimal\td__Archaea\t0\t0\t0.0\nminimal\td__Bacteria\t4\t0\t100.0\nanother\td__Archaea\t0\t0\t0.0\nanother\td__Bacteria\t0\t2\t0.0\ntotal\td__Archaea\t0\t0\t0.0\naverage\td__Archaea\t0.0\t0.0\tnan\ntotal\td__Bacteria\t4\t2\t66.7\naverage\td__Bacteria\t2.0\t1.0\t50.0\n", to_print.getvalue())

        to_print = StringIO()
        found_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print,
                                  binned_otu_table_io=found_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io)
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'])
                          ])+"\n",
                         found_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'])])+"\n",
                         not_found_otu_table_io.getvalue())

    def test_print_assembly_appraise(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        assembly_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        assemblies = "\n".join(["\t".join(x) for x in assembly_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        assembly_collection = OtuTableCollection()
        assembly_collection.add_otu_table(StringIO(assemblies))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 assembly_otu_table_collection=assembly_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)

        to_print = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print, doing_assembly=True)
        self.assertEqual("sample\tdomain\tnum_binned\tnum_assembled\tnum_not_found\tpercent_binned\tpercent_assembled\nminimal\td__Archaea\t0\t0\t0\t0.0\t0.0\nminimal\td__Bacteria\t0\t4\t0\t0.0\t100.0\nanother\td__Archaea\t0\t0\t0\t0.0\t0.0\nanother\td__Bacteria\t0\t0\t2\t0.0\t0.0\ntotal\td__Archaea\t0\t0\t0\t0.0\t0.0\naverage\td__Archaea\t0.0\t0.0\t0.0\tnan\tnan\ntotal\td__Bacteria\t0\t4\t2\t0.0\t66.7\naverage\td__Bacteria\t0.0\t2.0\t1.0\t0.0\t50.0\n", to_print.getvalue())

        to_print = StringIO()
        found_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print,
                                  doing_assembly=True,
                                  assembled_otu_table_io=found_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io)
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'])
                          ])+"\n",
                         found_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'])])+"\n",
                         not_found_otu_table_io.getvalue())

        # Check that unbinned is the same as assembled OTUs when no binning is done
        to_print = StringIO()
        assembled_otu_table_io = StringIO()
        unbinned_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print,
                                  doing_assembly=True,
                                  assembled_otu_table_io=assembled_otu_table_io,
                                  unbinned_otu_table_io=unbinned_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io)
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'])
                          ])+"\n",
                         assembled_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'])
                          ])+"\n",
                         unbinned_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'])])+"\n",
                         not_found_otu_table_io.getvalue())


    def test_print_assembly_appraise_all_binned(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        assembly_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        assemblies = "\n".join(["\t".join(x) for x in assembly_otu_table])

        genomes = assemblies

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        assembly_collection = OtuTableCollection()
        assembly_collection.add_otu_table(StringIO(assemblies))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 assembly_otu_table_collection=assembly_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)
        a = res[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_not_found)

        to_print = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print, doing_assembly=True)
        self.assertEqual("sample\tdomain\tnum_binned\tnum_assembled\tnum_not_found\tpercent_binned\tpercent_assembled\nminimal\td__Archaea\t0\t0\t0\t0.0\t0.0\nminimal\td__Bacteria\t4\t4\t0\t100.0\t100.0\nanother\td__Archaea\t0\t0\t0\t0.0\t0.0\nanother\td__Bacteria\t0\t0\t2\t0.0\t0.0\ntotal\td__Archaea\t0\t0\t0\t0.0\t0.0\naverage\td__Archaea\t0.0\t0.0\t0.0\tnan\tnan\ntotal\td__Bacteria\t4\t4\t2\t66.7\t66.7\naverage\td__Bacteria\t2.0\t2.0\t1.0\t50.0\t50.0\n", to_print.getvalue())

        # Check that unbinned is the same as assembled OTUs when no binning is done
        to_print = StringIO()
        assembled_otu_table_io = StringIO()
        unbinned_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print,
                                  doing_assembly=True,
                                  assembled_otu_table_io=assembled_otu_table_io,
                                  unbinned_otu_table_io=unbinned_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io)
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'])
                          ])+"\n",
                         assembled_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          ])+"\n",
                         unbinned_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers),
                          "\t".join(['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'])])+"\n",
                         not_found_otu_table_io.getvalue())

    def test_print_appraisal_output_found_in(self):
        metagenome_otu_table = [self.headers,
                    ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']
                    ]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])

        genomes_otu_table = [self.headers,['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']
                    ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 output_found_in=True,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[1]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual('genome', a.binned_otus[0].data[a.binned_otus[0].fields.index('found_in')])
        a = res[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual('', a.not_found_otus[0].data[a.not_found_otus[0].fields.index('found_in')])

        to_print = StringIO()
        found_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print,
                                  binned_otu_table_io=found_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io,
                                  output_found_in=True)
        self.assertEqual("\n".join([
                          "\t".join(self.headers + ['found_in']),
                          "\t".join(['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales', 'genome'])
                          ])+"\n",
                         found_otu_table_io.getvalue())
        self.assertEqual("\n".join([
                          "\t".join(self.headers + ['found_in']),
                          "\t".join(['4.11.ribosomal_protein_L10','another','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus', ''])])+"\n",
                         not_found_otu_table_io.getvalue())

    def test_print_appraisal_archive_input_output(self):
        version_4_fields = '"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?", "read_unaligned_sequences", "equal_best_hit_taxonomies", "taxonomy_assignment_method"]'
        alignment_shas = '"alignment_hmm_sha256s": ["4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9"]'
        package_shas = '"singlem_package_sha256s": ["e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48", "e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48", "e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48", "e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48"]'
        binned_otu = '["4.12.ribosomal_protein_L11_rplK", "minimal", "CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false, ["ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA", "ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA", "ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA", "ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA"], ["Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales"], "singlem_query_based"]'
        unbinned_otu = '["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false, ["ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA"], ["Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus"], "singlem_query_based"]'

        metagenomes = "".join([
            '{',
            ", ".join([
                '"version": 4',
                alignment_shas,
                package_shas,
                version_4_fields,
                "".join(['"otus": [', unbinned_otu, ',', binned_otu, ']']),
            ]),
            '}',
        ])

        genomes = "".join([
            '{',
            ", ".join([
                '"version": 4',
                alignment_shas,
                package_shas,
                version_4_fields,
                "".join(['"otus": [', binned_otu, ']']),
            ]),
            '}',
        ])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_archive_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_archive_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[0]
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(9.76/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(4.88/4)}, a.num_not_found)

        to_print = StringIO()
        found_otu_table_io = StringIO()
        not_found_otu_table_io = StringIO()
        appraiser.print_appraisal(app, packages, True, to_print,
                                  binned_otu_table_io=found_otu_table_io,
                                  unaccounted_for_otu_table_io=not_found_otu_table_io,
                                  output_style='archive'
                                  )

        expected_found = "".join([
            '{',
            ", ".join([
                '"version": 4',
                alignment_shas,
                package_shas,
                version_4_fields,
                "".join(['"otus": [', binned_otu, ']']),
            ]),
            '}',
        ])

        expected_not_found = "".join([
            '{',
            ", ".join([
                '"version": 4',
                alignment_shas,
                package_shas,
                version_4_fields,
                "".join(['"otus": [', unbinned_otu, ']']),
            ]),
            '}',
        ])

        self.assertEqual(expected_found, found_otu_table_io.getvalue())
        self.assertEqual(expected_not_found, not_found_otu_table_io.getvalue())

    def test_print_appraisal_archive_input_output_cli(self):
        with in_tempdir():
            cmd = (
                "{} appraise "
                "--metagenome-archive-otu-tables {}/appraise_example3/reads.otu_table.json "
                "--genome-otu-tables {}/appraise_example3/bins.otu_table.tsv "
                "--metapackage {}/four_package.smpkg "
                "--output-binned-otu-table binned.otu_table.json "
                "--output-unaccounted-for-otu-table unbinned.otu_table.json "
                "--output-style archive "
            ).format(
                path_to_script,
                path_to_data,
                path_to_data,
                path_to_data
            )
            extern.run(cmd)

            version_4_fields = '"fields": ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "read_names", "nucleotides_aligned", "taxonomy_by_known?", "read_unaligned_sequences", "equal_best_hit_taxonomies", "taxonomy_assignment_method"]'
            alignment_shas = '"alignment_hmm_sha256s": ["4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9", "4b0bf5b3d7fd2ca16e54eed59d3a07eab388f70f7078ac096bf415f1c04731d9"]'
            package_shas = '"singlem_package_sha256s": ["e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48", "e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48", "e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48", "e4de3077fe4f7869ae1d9c49fc650c664153325fd2bc5997044c983dedd36a48"]'
            binned_otu = '["4.12.ribosomal_protein_L11_rplK", "minimal", "GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC", 4, 9.75609756097561, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales", ["HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910", "HWI-ST1243:156:D1K83ACXX:7:1103:21187:63124", "HWI-ST1243:156:D1K83ACXX:7:1108:10813:6928", "HWI-ST1243:156:D1K83ACXX:7:1105:12385:81842"], [60, 60, 60, 60], false, ["ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA", "ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA", "ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA", "ATTAACAGTAGCTGAAGTTACTGACCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGCGTCGTGCAGCTGAA"], ["Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales"], "singlem_query_based"]'
            unbinned_otu = '["4.11.ribosomal_protein_L10", "minimal", "TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA", 2, 4.878048780487805, "Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus", ["HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482", "HWI-ST1243:156:D1K83ACXX:7:1105:19152:28331"], [60, 60], false, ["ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA"], ["Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus"], "singlem_query_based"]'
            expected_found = "".join([
                '{',
                ", ".join([
                    '"version": 4',
                    alignment_shas,
                    package_shas,
                    version_4_fields,
                    "".join(['"otus": [', binned_otu, ']']),
                ]),
                '}',
            ])

            expected_not_found = "".join([
                '{',
                ", ".join([
                    '"version": 4',
                    alignment_shas,
                    package_shas,
                    version_4_fields,
                    "".join(['"otus": [', unbinned_otu, ']']),
                ]),
                '}',
            ])

            with open('binned.otu_table.json') as f:
                observed_found = "".join(f.readlines())
            with open('unbinned.otu_table.json') as f:
                observed_not_found = "".join(f.readlines())

            self.assertEqual(expected_found, observed_found)
            self.assertEqual(expected_not_found, observed_not_found)

    def test_contamination(self):
        metagenome_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'minimal',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            [
                '4.16.ribosomal_protein_S5',
                'minimal',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '8','12.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        genomes_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC', #one base pair different to the one above
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(12.50/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_not_found)


    def test_contamination_near_enough(self):
        metagenome_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'minimal',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            [
                '4.16.ribosomal_protein_S5',
                'another',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '8','12.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        genomes_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC', #one base pair different to the one above
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGGGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',#one base pair different to the one above
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[0]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(12.50/4)}, a.num_not_found)
        a = res[1]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_not_found)

        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 sequence_identity=0.9,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        def compare_res(res): return res.metagenome_sample_name
        sorted_results = list(sorted(app.appraisal_results, key=compare_res))
        a = sorted_results[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(12.50/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)
        a = sorted_results[1]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_not_found)

    def test_for_concatenated_genomes(self):
        metagenome_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'minimal',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            [
                '4.16.ribosomal_protein_S5',
                'another',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '8','12.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        genomes_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGGGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',#one base pair different to the one above
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGGGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',#one base pair different to the one above
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = sorted(app.appraisal_results)
        a = res[0]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(12.50/4)}, a.num_not_found)
        a = res[1]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_not_found)

    def test_assembly_input(self):
        metagenome_otu_table = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','minimal','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus'],
            ['4.14.ribosomal_protein_L16_L10E_rplP','minimal','CAAAAAAAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','5','10.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        assembly_otu_table = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','assembly','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.007','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.ribosomal_protein_L10','assembly','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','1.01','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus']]
        assemblies = "\n".join(["\t".join(x) for x in assembly_otu_table])

        genomes_otu_table = [
            self.headers,
            ['4.12.ribosomal_protein_L11_rplK','genome','GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC','1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli']]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        assembly_collection = OtuTableCollection()
        assembly_collection.add_otu_table(StringIO(assemblies))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 assembly_otu_table_collection=assembly_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(1, len(app.appraisal_results))
        a = app.appraisal_results[0]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(26.83/4)}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(10.76/4)}, a.num_not_found)
        self.assertEqual('minimal', a.metagenome_sample_name)
        self.assertEqual(1, len(a.binned_otus))
        self.assertEqual('GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                         a.binned_otus[0].sequence)
        self.assertEqual(2, len(a.assembled_otus))
        self.assertEqual(sorted(['GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                                 'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG']),
                         sorted([o.sequence for o in a.assembled_otus]))
        self.assertEqual(1, len(a.not_found_otus))
        self.assertEqual('CAAAAAAAAAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                         a.not_found_otus[0].sequence)

    def test_appraise_plot_real_data(self):
        """Not a real test, just developing the code"""
        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        with open(os.path.join(path_to_data, 'appraise_example4', 'SRR5040536.reads.long_sample_names.otu_table.csv')) as f:
            metagenome_collection.add_otu_table(f)
        genome_collection = OtuTableCollection()
        with open(os.path.join(path_to_data, 'appraise_example4', 'SRR5040536.binned.otu_table.csv')) as f:
            genome_collection.add_otu_table(f)
        assembly_collection = OtuTableCollection()
        with open(os.path.join(path_to_data, 'appraise_example4', 'SRR5040536.assembly.otu_table.csv')) as f:
            assembly_collection.add_otu_table(f)
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 assembly_otu_table_collection=assembly_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)

        with tempfile.NamedTemporaryFile(mode='w',suffix='.svg',prefix='single_test_appraisal.') as f:
            app.plot(
                output_svg_base=f.name,
                cluster_identity = 0.89,
                doing_assembly=True,
                doing_binning=True
            )

    def test_appraise_assembly_imperfectly(self):
        metagenome_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'minimal',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '7','17.07','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            [
                '4.16.ribosomal_protein_S5',
                'another',
                'GGTACCGGCGTCATCGCCGGTGGCGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '8','12.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria'],
            [
                '4.16.ribosomal_protein_S5',
                'another',
                'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGGCCGGCATCCGC', # way different to the one above
                '9','17.50','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']]
        metagenomes = "\n".join(["\t".join(x) for x in metagenome_otu_table])
        genomes_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'genome',
                'AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC', #one base pair different to the one above
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'genome',
                'GGTACCGGCGTCATCGCCGGTGGGGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        genomes = "\n".join(["\t".join(x) for x in genomes_otu_table])
        assemblies_otu_table = [
            self.headers,[
                '4.12.ribosomal_protein_L11_rplK',
                'assembly',
                'GGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC',
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.12.ribosomal_protein_L11_rplK',
                'assembly',
                'AGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATC', #one base pair different to the one above
                '1','1.02','Root; d__Bacteria; p__Firmicutes; c__Bacilli'],
            [
                '4.16.ribosomal_protein_S5',
                'assembly',
                'GGTACCGGCGTCATCGCCGGTGGGGCGGCACGCGCCATCTTGGAGATGGCCGGCATCCGC',
                '1','1.06','Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria']
        ]
        assemblies = "\n".join(["\t".join(x) for x in assemblies_otu_table])

        appraiser = Appraiser()
        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table(StringIO(metagenomes))
        genome_collection = OtuTableCollection()
        genome_collection.add_otu_table(StringIO(genomes))
        assembly_collection = OtuTableCollection()
        assembly_collection.add_otu_table(StringIO(assemblies))
        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 assembly_otu_table_collection=assembly_collection,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = self._sort_appraisal_results(app.appraisal_results)
        a = res[0]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round((12.5+17.5)/4)}, a.num_not_found)
        a = res[1]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)

        packages = Metapackage.acquire(os.path.join(path_to_data, 'four_package.smpkg')).singlem_packages
        app = appraiser.appraise(genome_otu_table_collection=genome_collection,
                                 metagenome_otu_table_collection=metagenome_collection,
                                 assembly_otu_table_collection=assembly_collection,
                                 sequence_identity=0.9,
                                 packages=packages,
                                 window_size=DEFAULT_WINDOW_SIZE)
        self.assertEqual(2, len(app.appraisal_results))
        res = self._sort_appraisal_results(app.appraisal_results)
        a = res[0]
        self.assertEqual('another', a.metagenome_sample_name)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(12.5/4)}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(12.5/4)}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.5/4)}, a.num_not_found)
        a = res[1]
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_binned)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': round(17.07/4)}, a.num_assembled)
        self.assertEqual({'d__Archaea': 0, 'd__Bacteria': 0}, a.num_not_found)


if __name__ == "__main__":
    unittest.main()
