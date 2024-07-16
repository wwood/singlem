
#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft, Samuel Aroney.
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
import os
import sys
from bird_tool_utils import in_tempdir
import extern

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.regenerator import Regenerator
from singlem.singlem_package import SingleMPackage

# small package with few sequences
input_package = os.path.join(path_to_data, "4.11.22seqs.v3.gpkg.spkg")
# same sequences with altered taxonomy (rotated by 1)
input_sequences = os.path.join(path_to_data, "regenerate", "input_sequences.fasta")
input_taxonomy = os.path.join(path_to_data, "regenerate", "input_taxonomy.tsv")
# few sequences from uniprot that match input package
euk_sequences = os.path.join(path_to_data, "regenerate", "uniprot_sprot_hits.fa")
euk_taxonomy = os.path.join(path_to_data, "regenerate", "uniprot_sprot_taxonomy.tsv")

output_package = "S3.regenerated.gpkg.spkg"
sequence_prefix = "prefix~"
window_position = 20
min_aligned_percent = 30

class Tests(unittest.TestCase):
    def test_regenerate_package(self):
        self.maxDiff = None
        with in_tempdir():
            Regenerator().regenerate(
                    input_singlem_package = input_package,
                    input_sequences = input_sequences,
                    input_taxonomy = input_taxonomy,
                    euk_sequences = euk_sequences,
                    euk_taxonomy = euk_taxonomy,
                    window_position = window_position,
                    output_singlem_package = output_package,
                    sequence_prefix = sequence_prefix,
                    min_aligned_percent = min_aligned_percent,
                    no_further_euks = False)

            pkg = SingleMPackage.acquire(output_package)
            self.assertEqual(window_position, pkg.singlem_position())

            # assert sequences and taxonomy have been supplemented with euk sequences, updated and trimmed
            # observed_output_fasta = list(io.open(pkg.graftm_package().unaligned_sequence_database_path()))
            with open(pkg.graftm_package().unaligned_sequence_database_path()) as f:
                observed_output_fasta = list(f)
            # expected_output_fasta = list(io.open(os.path.join(path_to_data, "regenerate", "output_trimmed.fasta")))
            with open(os.path.join(path_to_data, "regenerate", "output_trimmed.fasta")) as f:
                expected_output_fasta = list(f)
            self.assertListEqual(observed_output_fasta, expected_output_fasta)

            # observed_output_seqinfo = list(io.open(pkg.graftm_package().taxtastic_seqinfo_path()))
            with open(pkg.graftm_package().taxtastic_seqinfo_path()) as f:
                observed_output_seqinfo = list(f)
            # expected_output_seqinfo = list(io.open(os.path.join(path_to_data, "regenerate", "output_seqinfo.csv")))
            with open(os.path.join(path_to_data, "regenerate", "output_seqinfo.csv")) as f:
                expected_output_seqinfo = list(f)
            self.assertListEqual(sorted(observed_output_seqinfo), sorted(expected_output_seqinfo))

            # assert seed_idx file has been created
            seed_idx = pkg.graftm_package().diamond_database_path() + ".seed_idx"
            self.assertTrue(os.path.exists(seed_idx))

            # assert taxonomy hash file has been created, including duplicate seq removed by graftM
            observed_taxonomy_hash = pkg.taxonomy_hash()
            expected_taxonomy_hash = {
            'prefix~1': ['d__Bacteria', 'p__Cyanobacteria', 'c__[Chroococcales]', 'o__Chroococcales', 'f__[Synechococcus]', 'g__Synechococcus', 's__Synechococcus_sp.'],
            'prefix~2': ['d__Bacteria', 'p__Proteobacteria', 'c__Deltaproteobacteria', 'o__Myxococcales', 'f__Myxococcaceae', 'g__Myxococcus', 's__Myxococcus_stipitatus'],
            'prefix~3': ['d__Bacteria', 'p__Proteobacteria', 'c__Betaproteobacteria', 'o__Nitrosomonadales', 'f__Nitrosomonadaceae', 'g__Nitrosomonas', 's__Nitrosomonas_cryotolerans'],
            'prefix~4': ['d__Archaea', 'p__Crenarchaeota', 'c__Thermoprotei', 'o__Desulfurococcales', 'f__Desulfurococcaceae', 'g__Thermosphaera', 's__Thermosphaera_aggregans'],
            'prefix~5': ['d__Bacteria', 'p__Proteobacteria', 'c__Betaproteobacteria', 'o__Rhodocyclales', 'f__Rhodocyclaceae', 'g__Azoarcus', 's__Azoarcus_toluclasticus'],
            'prefix~6': ['d__Bacteria', 'p__Proteobacteria', 'c__Betaproteobacteria', 'o__[Candidatus_Accumulibacter]', 'f__[Candidatus_Accumulibacter]', 'g__Candidatus_Accumulibacter', 's__Candidatus_Accumulibacter_sp._BA-92'],
            'prefix~7': ['d__Bacteria', 'p__Cyanobacteria', 'c__[Chroococcales]', 'o__Chroococcales', 'f__[Acaryochloris]', 'g__Acaryochloris', 's__Acaryochloris_sp._CCMEE_5410'],
            'prefix~8': ['d__Bacteria', 'p__Bacteroidetes', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Prevotellaceae', 'g__Prevotella', 's__Prevotella_corporis'],
            'prefix~9': ['d__Bacteria', 'p__Cyanobacteria', 'c__[Prochlorales]', 'o__Prochlorales', 'f__Prochloraceae', 'g__Prochloron', 's__Prochloron_didemni'],
            'prefix~10': ['d__Bacteria', 'p__Bacteroidetes', 'c__Cytophagia', 'o__Cytophagales', 'f__Cyclobacteriaceae', 'g__Indibacter', 's__Indibacter_alkaliphilus'],
            'prefix~11': ['d__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Sphingomonadales', 'f__Sphingomonadaceae', 'g__Sphingobium', 's__Sphingobium_sp._KK22'],
            'prefix~12': ['d__Bacteria', 'p__Firmicutes', 'c__Clostridia', 'o__Clostridiales', 'f__Lachnospiraceae', 'g__Blautia', 's__Blautia_wexlerae'],
            'prefix~13': ['d__Bacteria', 'p__Firmicutes', 'c__Clostridia', 'o__Clostridiales', 'f__Lachnospiraceae', 'g__[Lachnospiraceae_bacterium_NK4A179]', 's__Lachnospiraceae_bacterium_NK4A179'],
            'prefix~14': ['d__Bacteria', 'p__Actinobacteria', 'c__Actinobacteria', 'o__Actinomycetales', 'f__Microbacteriaceae', 'g__Microbacterium', 's__Microbacterium_sp.'],
            'prefix~15': ['d__Bacteria', 'p__Firmicutes', 'c__Bacilli', 'o__Lactobacillales', 'f__Carnobacteriaceae', 'g__Carnobacterium', 's__Carnobacterium_gallinarum'],
            'prefix~16': ['d__Bacteria', 'p__Proteobacteria', 'c__Alphaproteobacteria', 'o__Rhodospirillales', 'f__Acetobacteraceae', 'g__Kozakia', 's__Kozakia_baliensis'],
            'prefix~17': ['d__Archaea', 'p__Crenarchaeota', 'c__Thermoprotei', 'o__Thermoproteales', 'f__Thermoproteaceae', 'g__Pyrobaculum', 's__Pyrobaculum_calidifontis'],
            'prefix~18': ['d__Bacteria', 'p__Cyanobacteria', 'c__[Nostocales]', 'o__Nostocales', 'f__Nostocaceae', 'g__Anabaena', 's__Anabaena_sp._90'],
            'prefix~19': ['d__Bacteria', 'p__Bacteroidetes', 'c__Flavobacteriia', 'o__Flavobacteriales', 'f__Flavobacteriaceae', 'g__Aquimarina', 's__Aquimarina_sp._SW150'],
            'prefix~20': ['d__Archaea', 'p__Crenarchaeota', 'c__Thermoprotei', 'o__Thermoproteales', 'f__Thermoproteaceae', 'g__Vulcanisaeta', 's__Vulcanisaeta_moutnovskia'],
            'prefix~21': ['d__Archaea', 'p__Crenarchaeota', 'c__Thermoprotei', 'o__Thermoproteales', 'f__Thermoproteaceae', 'g__Vulcanisaeta', 's__Vulcanisaeta_moutnovskia'],
            'prefix~RK10_TOBAC': ['d__Eukaryota', 'Viridiplantae'],
            'prefix~RK10_ARATH': ['d__Eukaryota', 'Viridiplantae'],
            'prefix~RK10_SPIOL': ['d__Eukaryota', 'Viridiplantae'],
            'prefix~RM11_YEAST': ['d__Eukaryota', 'Fungi'],
            'prefix~RM10_DROME': ['d__Eukaryota', 'Metazoa'],
            'prefix~RM10_MOUSE': ['d__Eukaryota', 'Metazoa'],
            'prefix~RM10_RAT': ['d__Eukaryota', 'Metazoa'],
            'prefix~RM10_DROPS': ['d__Eukaryota', 'Metazoa'],
            'prefix~RM10_DANRE': ['d__Eukaryota', 'Metazoa'],
            'prefix~RM10_BOVIN': ['d__Eukaryota', 'Metazoa'],
            'prefix~RM10_HUMAN': ['d__Eukaryota', 'Metazoa'],
            'prefix~RLA0_YEAST': ['d__Eukaryota', 'Fungi'],
            'prefix~RLA0_PODAS': ['d__Eukaryota', 'Fungi'],
            'prefix~RLA0_SCHPO': ['d__Eukaryota', 'Fungi'],
            }
            self.assertEqual(observed_taxonomy_hash, expected_taxonomy_hash)


    #@unittest.skip("CLI testing is so slow. Can't figure out how to mock with extern.")
    def test_hello_word_cmdline(self):
        with in_tempdir():
            cmd = "{} regenerate ".format(path_to_script)
            cmd += "--input-singlem-package {} ".format(input_package)
            cmd += "--input-sequences {} ".format(input_sequences)
            cmd += "--input-taxonomy {} ".format(input_taxonomy)
            cmd += "--candidate-decoy-sequences {} ".format(euk_sequences)
            cmd += "--candidate-decoy-taxonomy {} ".format(euk_taxonomy)
            cmd += "--window-position {} ".format(window_position)
            cmd += "--sequence-prefix {} ".format(sequence_prefix)
            cmd += "--output-singlem-package {} ".format(output_package)
            extern.run(cmd)

            SingleMPackage.acquire(output_package)


if __name__ == "__main__":
    unittest.main()
