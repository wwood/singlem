#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2023 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2023"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os

import extern
from bird_tool_utils import in_tempdir

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from singlem.metapackage import Metapackage, CUSTOM_TAXONOMY_DATABASE_NAME

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument(
        '--metapackage1',
        help="First input smpkg", required=True)
    parent_parser.add_argument(
        '--metapackage2',
        help="Second input smpkg", required=True)
    parent_parser.add_argument(
        '--output-metapackage',
        help="Output smpkg", required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    smpkg1 = Metapackage.acquire(os.path.abspath(args.metapackage1))
    smpkg2 = Metapackage.acquire(os.path.abspath(args.metapackage2))
    output_metapackage = os.path.abspath(args.output_metapackage)

    # in tempdir
    with in_tempdir():
        # Create a new fasta file containing the prefilter seqs from metapackage
        logging.info("Extracting prefilter seqs ..")
        prefilter1_seqs = 'prefilter1_seqs.fasta'
        extern.run('diamond getseq -d {0} -o {1}'.format(smpkg1.prefilter_db_path(), prefilter1_seqs))
        prefilter2_seqs = 'prefilter2_seqs.fasta'
        extern.run('diamond getseq -d {0} -o {1}'.format(smpkg2.prefilter_db_path(), prefilter2_seqs))

        # Create diamond db from prefilter seqs
        logging.info("Creating diamond db ..")
        extern.run(f"cat {prefilter1_seqs} {prefilter2_seqs} > prefilter_seqs.fasta")
        extern.run(f"diamond makedb --in prefilter_seqs.fasta -d prefilter_seqs.dmnd")

        logging.info("Generating metapackage ..")
        logging.warning("Merged package does not contain nucleotide SDB or taxon genome lengths")
        Metapackage.generate(
            singlem_packages=[s.base_directory() for s in (smpkg1.singlem_packages + smpkg2.singlem_packages)],
            nucleotide_sdb=None,
            prefilter_clustering_threshold=None,
            output_path=output_metapackage,
            threads=1,
            prefilter_diamond_db='prefilter_seqs.dmnd',
            taxon_genome_lengths=None,
            taxonomy_database_name=f'metapackage merged from {args.metapackage1} and {args.metapackage2}',
            taxonomy_database_version=None,
            diamond_prefilter_performance_parameters=None,
            diamond_taxonomy_assignment_performance_parameters=None,
            makeidx_sensitivity_params=None,
        )

        logging.info("Done")
