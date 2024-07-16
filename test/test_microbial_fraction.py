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
import extern

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path


class Tests(unittest.TestCase):
    output_headers = str.split('sample  bacterial_archaeal_bases        metagenome_size   read_fraction    average_bacterial_archaeal_genome_size  warning')
    maxDiff = None

    def test_marine0(self):
        cmd = "{} microbial_fraction -p {}/read_fraction/marine0.profile  --input-metagenome-sizes {}/read_fraction/marine0.num_bases --taxon-genome-lengths-file {}/read_fraction/gtdb_mean_genome_sizes.tsv".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('\t'.join(self.output_headers)+'\n' + '\t'.join(str.split('marine0.1       16593586562      17858646300.0   92.92   3718692')+[''])+'\n', obs)

    def test_smafa_count_unpaired(self):
        cmd = "{} microbial_fraction -p {}/read_fraction/marine0.profile  --forward {}/read_fraction/marine0.1.fa --taxon-genome-lengths-file {}/read_fraction/gtdb_mean_genome_sizes.tsv".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('\t'.join(self.output_headers)+'\n' + '\t'.join(str.split('marine0.1       16593586562      3       100.00  3718692')+['WARNING: The most abundant taxons not assigned to the species level account for a large fraction of the total estimated read fraction. This may mean that the read_fraction estimate is inaccurate.'])+'\n', obs)

    def test_smafa_count_paired(self):
        cmd = "{} microbial_fraction -p {}/read_fraction/marine0.profile  --forward {}/read_fraction/marine0.1.fa --reverse {}/read_fraction/marine0.2.fa --taxon-genome-lengths-file {}/read_fraction/gtdb_mean_genome_sizes.tsv".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('\t'.join(self.output_headers)+'\n' + '\t'.join(str.split('marine0.1       16593586562      6       100.00    3718692')+['WARNING: The most abundant taxons not assigned to the species level account for a large fraction of the total estimated read fraction. This may mean that the read_fraction estimate is inaccurate.'])+'\n', obs)

    def test_output_per_taxon_read_fractions(self):
        cmd = "{} microbial_fraction -p <(head -5 {}/read_fraction/marine0.profile) --input-metagenome-sizes {}/read_fraction/marine0.num_bases --taxon-genome-lengths-file {}/read_fraction/gtdb_mean_genome_sizes.tsv --output-tsv /dev/null --output-per-taxon-read-fractions /dev/stdout".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('''sample	taxonomy	base_contribution
marine0.1	d__Archaea	5608320.893634946
marine0.1	p__Thermoproteota	858552.2612056588
marine0.1	p__Desulfobacterota	2716953.844763834
marine0.1	p__Proteobacteria	7151244.856821437
''', obs)

    def test_average_genome_size(self):
        '''Here the ave genome size has been checked by hand'''

        cmd = "{} microbial_fraction -p <(head -5 {}/read_fraction/marine0.profile)  --input-metagenome-sizes {}/read_fraction/marine0.num_bases --taxon-genome-lengths-file {}/read_fraction/gtdb_mean_genome_sizes.tsv".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data)
        obs = extern.run(cmd)
        self.assertEqual('\t'.join(self.output_headers)+'\n' + '\t'.join(str.split('marine0.1       16335072      17858646300.0   0.09  2278253')+[''])+'\n', obs)




if __name__ == "__main__":
    unittest.main()
