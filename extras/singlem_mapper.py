#!/usr/bin/env python3.6

###############################################################################
#
#    Copyright (C) 2020 Ben Woodcroft
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
__copyright__ = "Copyright 2022"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import os
import sys
import argparse
import logging
import tempfile

import extern
import datatable as dt
from datatable import join, f

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--forward', help='forward reads for mapping', required=True)
    parent_parser.add_argument('--reverse', help='reverse reads for mapping', required=True)
    parent_parser.add_argument('--singlem-profile', help='regular or output extras form of singlem profile for the reads being mapped', required=True)
    parent_parser.add_argument('--reference-genome-singlem-database', help='map reads only to those genomes closely related in this SingleM database', required=True)
    parent_parser.add_argument('--reference-genome-locations', help='Two column TSV file with <genome_id>TAB<location>', required=True)
    parent_parser.add_argument('--output', help='output profile from CoverM', required=True)
    parent_parser.add_argument('--coverm-extra-args', help='extra arguments for CoverM e.g. --threads')

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read reference genome locations
    logging.info("Reading reference genome locations")
    genome_locations = dt.fread(args.reference_genome_locations, header=False, columns=['sample', 'location'])
    genome_locations[:, 'location'] = [os.path.join(os.path.dirname(args.reference_genome_locations), x) for x in genome_locations['location'].to_list()[0]]

    # Query profile against genome DB
    logging.info("Querying genome database with SingleM query ..")
    cmd = "singlem query --query-otu-table {} --db {} --max-divergence {}".format(
        args.singlem_profile, args.reference_genome_singlem_database, 0)
    query_result = dt.fread(extern.run(cmd))

    # Find reference genome locations that match the read profile
    hit_samples = dt.unique(query_result[:,f.sample])
    hit_samples[:,'present']=True
    hit_samples.key = 'sample'
    g2 = genome_locations[:,:,join(hit_samples)][f.present==1,:][:,[f.location]]
    if g2.nrows != hit_samples.nrows:
        raise Exception("Some samples that are present the sample DB were not found in the genome database")
        
    logging.info("Found %d matching genomes" % g2.nrows)
    # Run CoverM 
    logging.info("Running CoverM ..")
    with tempfile.NamedTemporaryFile(mode='w', prefix='singlem_mapper_genome_list') as f:
        g2.to_csv(f.name, header=False)
        f.flush()

        cmd = "coverm genome --genome-fasta-list {} -1 {} -2 {} -o {} {}".format(
            f.name, args.forward, args.reverse, args.output, args.coverm_extra_args if args.coverm_extra_args else '')
        extern.run(cmd)
    logging.info("Finished")
