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
import shutil

import extern
from bird_tool_utils import in_tempdir

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from singlem.metapackage import Metapackage

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument(
        '--metapackage',
        help="Input smpkg", required=True)
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

    logging.info("Loading metapackage from %s", args.metapackage)
    smpkg = Metapackage.acquire(os.path.abspath(args.metapackage))
    output_metapackage = os.path.abspath(args.output_metapackage)

    # Copy input metapackage to output
    logging.info("Copying metapackage to %s", output_metapackage)
    shutil.copytree(args.metapackage, output_metapackage)

    logging.info("Loading metapackage from %s", output_metapackage)
    metapackage = Metapackage.acquire(output_metapackage)

    # Delete read_taxonomies
    logging.info("Deleting various files")
    os.remove(metapackage._sqlite_db_path)

    if metapackage.nucleotide_sdb_path() is not None:
        logging.info("Deleting nucleotide sdb")
        shutil.rmtree(metapackage.nucleotide_sdb_path())

    # From each graftm package
    for singlem_package in metapackage.singlem_packages:
        os.remove(singlem_package.taxonomy_hash_path())

        gpkg = singlem_package.graftm_package()

        # Delete unaligned fasta
        shutil.rmtree(gpkg.reference_package_path())

        os.remove(gpkg.diamond_database_path())
        os.remove(gpkg.diamond_database_path() + ".seed_idx")

        # Bit of a hack - the seq IDs are read from the unaligned fasta file,
        # but the sequences are not used. So here we just make a dummy fasta
        # file with 1 char sequences to save space.
        seq_ids = singlem_package.get_sequence_ids()
        with open(gpkg.unaligned_sequence_database_path(), "w") as f:
            for seq_id in seq_ids:
                f.write(">%s\nA\n" % seq_id)
