
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

from singlem.singlem_package import SingleMPackage
import unittest
from unittest.mock import Mock,patch
import os
import sys
from bird_tool_utils import in_tempdir
import extern
import io

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.regenerator import Regenerator

class Tests(unittest.TestCase): 
    def test_regenerate_package(self):
        with in_tempdir():
            # small package with few sequences
            input_package = os.path.join(path_to_data, "4.11.22seqs.v3.gpkg.spkg")
            # same sequences with altered taxonomy (rotated by 1)
            input_sequences = os.path.join(path_to_data, "regenerate", "input_sequences.fasta")
            input_taxonomy = os.path.join(path_to_data, "regenerate", "input_taxonomy.tsv")
            # few sequences from uniprot that match input package
            euk_sequences = os.path.join(path_to_data, "regenerate", "uniprot_sprot_hits.fa")
            euk_taxonomy = os.path.join(path_to_data, "regenerate", "uniprot_sprot_taxonomy.tsv")

            output_package = "regenerated.spkg"
            window_position = 20
            min_aligned_percent = 30

            
            Regenerator().regenerate(
                    input_singlem_package = input_package,
                    input_sequences = input_sequences,
                    input_taxonomy = input_taxonomy,
                    euk_sequences = euk_sequences,
                    euk_taxonomy = euk_taxonomy,
                    output_singlem_package = output_package,
                    min_aligned_percent = min_aligned_percent)

            pkg = SingleMPackage.acquire(output_package)
            #self.assertEqual(window_position, pkg.singlem_position())

            # assert sequences and taxonomy have been supplemented with euk sequences and updated
            observed_output_fasta = list(io.open(os.path.join(output_package, "regenerated", "4.11.22seqs_final_sequences.faa")))
            expected_output_fasta = list(io.open(os.path.join(path_to_data, "regenerate", "output_full.fasta")))
            self.assertListEqual(observed_output_fasta, expected_output_fasta)

            observed_output_seqinfo = list(io.open(os.path.join(output_package, "regenerated", "4.11.22seqs_final.gpkg.refpkg", "4_seqinfo.csv")))
            expected_output_seqinfo = list(io.open(os.path.join(path_to_data, "regenerate", "output_seqinfo.csv")))
            self.assertListEqual(observed_output_seqinfo, expected_output_seqinfo)
    

    def test_hello_word_cmdline(self):
        with in_tempdir():
            output_package = "regenerated.spkg"

            cmd = "{} regenerate ".format(path_to_script)
            cmd += "--input_singlem_package {} ".format(os.path.join(path_to_data, "4.11.22seqs.v3.gpkg.spkg"))
            cmd += "--input_sequences {} ".format(os.path.join(path_to_data, "regenerate", "input_sequences.fasta"))
            cmd += "--input_taxonomy {} ".format(os.path.join(path_to_data, "regenerate", "input_taxonomy.tsv"))
            cmd += "--euk_sequences {} ".format(os.path.join(path_to_data, "regenerate", "uniprot_sprot_hits.fa"))
            cmd += "--euk_taxonomy {} ".format(os.path.join(path_to_data, "regenerate", "uniprot_sprot_taxonomy.tsv"))
            cmd += "--output_singlem_package {} ".format(output_package)
            cmd += "--min_aligned_percent 30 "
            extern.run(cmd)

            SingleMPackage.acquire(output_package)


if __name__ == "__main__":
    unittest.main()
