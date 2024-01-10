#!/usr/bin/env python3

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

import argparse
import logging
import sys
import os

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--sra-identifier',required=True)
    parser.add_argument('--metapackage',required=True)

    parser.add_argument('--sra-download-methods', default='aws-http prefetch --guess-aws-location')

    args = parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    def run_singlem(sequence_input_arg, output_file_basename, use_sra):
        cmd = f'singlem pipe {sequence_input_arg} --no-assign-taxonomy --archive-otu-table >(gzip >{output_file_basename}.json.gz) --threads 1 --metapackage {args.metapackage}'
        if args.debug:
            cmd += ' --debug'

        try:
            logging.info(f"Attempting SingleM command: {cmd}")
            extern.run(cmd)
            return True
        except extern.ExternCalledProcessError as e:
            if use_sra:
                logging.warning(e)
                return False
            else:
                logging.error(e)
                raise(e)

    def analyse(run, use_sra):
        if use_sra:
            sequence_input_arg = f'--sra-files {run}.sra'
            return run_singlem(sequence_input_arg, run, use_sra)
        else:
            some_ena_analysed = False
            paired_result = None
            unpaired_result = None
            if os.path.exists(f'{run}_1.fastq.gz'):
                if os.path.exists(f'{run}_2.fastq.gz'):
                    sequence_input_arg = f'--forward {run}_1.fastq.gz --reverse {run}_2.fastq.gz'
                    some_ena_analysed = True
                    paired_result = run_singlem(sequence_input_arg, f'{run}_paired', use_sra)
                else:
                    raise Exception("Found a forward read file but no reverse read file from ENA")
            
            if os.path.exists(f'{run}.fastq.gz'):
                sequence_input_arg = f'--forward {run}.fastq.gz'
                some_ena_analysed = True
                unpaired_result = run_singlem(sequence_input_arg, f'{run}_unpaired', use_sra)

            if not some_ena_analysed:
                raise Exception("Unexpected form of ENA download")

            if paired_result is None:
                return unpaired_result
            elif unpaired_result is None:
                return paired_result
            else:
                return paired_result and unpaired_result



    logging.info(f"Attempting download of SRA format {args.sra_identifier} ..")
    sra_download_worked = False
    ena_download_worked = False
    run = args.sra_identifier
    try:
        extern.run(f'kingfisher get -r {run} --output-format-possibilities sra --hide-download-progress -m {args.sra_download_methods}')
        logging.info("Kingfisher get of .sra format worked")
        sra_download_worked = analyse(run,True) # Analyse too in case download works but singlem fails
    except extern.ExternCalledProcessError as e:
        logging.warning(e)
        pass

    try:
        if not sra_download_worked:
            logging.info("Failed .sra format download, trying ENA ..")
            extern.run(f'kingfisher get -r {run} --output-format-possibilities fastq.gz --hide-download-progress -m ena-ftp')
            logging.info("Kingfisher get of FASTQ.GZ from ENA format worked")
            analyse(run,False)
            ena_download_worked = True
    except extern.ExternCalledProcessError as e:
        logging.error(e)
        pass

    # Delete files for cleanliness
    for f in [f'{run}.sra',f'{run}_1.fastq.gz',f'{run}_2.fastq.gz',f'{run}.fastq.gz']:
        try:
            if os.path.exists(f):
                os.remove(f)
        except:
            logging.debug(f'File {f} not found when deleting')
            pass

    if sra_download_worked:
        logging.info("SRA-method worked")
    elif ena_download_worked:
        logging.info("ENA-method worked")
    else:
        raise Exception("No method worked")

