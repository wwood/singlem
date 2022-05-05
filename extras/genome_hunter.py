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

import argparse
from enum import unique
import logging
import sys
import os
import pandas as pd
from io import StringIO
import re
import tempfile

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def find_contigs(args):
    # Read query OTU table
    # Remove any OTUs that are not in the sample list
    queries = pd.read_csv(args.query_against_reads,sep='\t')
    if args.sample_list is not None:
        with open(args.sample_list) as f:
            samples_of_interest = set([s.strip() for s in f.readlines()])
        logging.info("Found {} samples of interest".format(len(samples_of_interest)))
        target_sequences = set(pd.merge(pd.DataFrame({'sample':list(samples_of_interest)}), queries, on='sample')['hit_sequence'])
        logging.info("Found {} target sequences from those samples of interest".format(len(target_sequences)))
        if len(target_sequences) == 0:
            raise Exception("No target sequences found")
    else:
        target_sequences = set(queries['hit_sequence'])
    

    # Run SingleM pipe on assembly files, putting results in STDOUT and reading, I guess
    assembly_file_otu_tables = {}
    for assembly_file in args.assembly_files:
        logging.info("Running SingleM on {}".format(assembly_file))
        cmd = f'singlem pipe --threads {args.threads} --genome-fasta-files {assembly_file} --otu_table /dev/stdout --no-assign-taxonomy --output-extras --singlem-metapackage {args.metapackage}'
        logging.debug("Running {}".format(cmd))
        output = extern.run(cmd)
        assembly_file_otu_tables[assembly_file] = pd.read_csv(StringIO(output),sep='\t')
    logging.info("Finished running SingleM on {} assemblies".format(len(assembly_file_otu_tables)))

    # Find OTUs in assembly singlem that are in the trimmed query OTU table
    good_otus = None
    for assembly_file,otu_table in assembly_file_otu_tables.items():
        logging.debug("Inspecting OTU table for {}".format(assembly_file))
        found_otus = pd.merge(otu_table, pd.DataFrame({'sequence':list(target_sequences)}), on='sequence')
        logging.info("Found {} assembled matching OTUs in {}".format(len(found_otus), assembly_file))
        new_frame = pd.DataFrame(found_otus)
        new_frame['assembly_file'] = assembly_file

        if good_otus is None:
            good_otus = new_frame
        else:
            good_otus = good_otus.append(new_frame)

    logging.info("Found {} matching OTUs".format(len(good_otus)))

    # Print the length of the contigs that match OTU sequences in the query OTU table
    last_assembly_file = None
    last_assembly_file_contig_lengths = None
    r = re.compile('(.*)_\d+_\d+_\d+')
    for otu_tuple in good_otus.iterrows():
        otu = otu_tuple[1]

        if last_assembly_file is None or last_assembly_file != otu['assembly_file']:
            last_assembly_file = otu['assembly_file']
            logging.info("Finding contig lengths from {} ..".format(last_assembly_file))
            last_assembly_file_contig_lengths_str = extern.run("bbbin describe <{}".format(last_assembly_file))
            last_assembly_file_contig_lengths = pd.read_csv(StringIO(last_assembly_file_contig_lengths_str),sep='\t',header=None)

        # Read sequence length of each contig in the original assembly file - unaligned sequences in the OTU table are from the transcripts, not the contigs
        assembly_file = otu['assembly_file']

        sequence_names = otu['read_names'].split(' ')
        for name in sequence_names:
            contig_name = r.match(name).group(1)
            sequence_length = last_assembly_file_contig_lengths[last_assembly_file_contig_lengths[0] == contig_name][1].values[0]

            print("\t".join((otu['assembly_file'], str(sequence_length), '\t'.join([str(s) for s in otu.values]))))

def classify_contigs(args):
    logging.info("Reading GTDB metadata ..")
    gtdb_bacteria = pd.read_csv(args.gtdb_bacteria_metadata, sep="\t")
    gtdb_archaea = pd.read_csv(args.gtdb_archaea_metadata, sep="\t")
    gtdb_metadata = pd.concat([gtdb_bacteria, gtdb_archaea])


    with open(args.output_best_hits, 'w') as best_hits_output_f:
        find_contigs_df = pd.read_csv(args.find_contigs_output, sep='\t', header=None)
        passing_rows = find_contigs_df[find_contigs_df[1] >= args.min_contig_length]
        logging.info("Found {} contigs of length >= {}".format(len(passing_rows), args.min_contig_length))
        for (assembly, contig_hits) in passing_rows.groupby(0):
            logging.info("Processing {} ..".format(assembly))

            singlem_orf_names = contig_hits[8].values
            # Extract contigs into tempfile
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fasta') as contigs_f:

                with tempfile.NamedTemporaryFile(mode='w',suffix='.names') as contig_names_f:
                    contig_names = set()
                    r = re.compile('(.*)_\d+_\d+_\d+')
                    for names in singlem_orf_names:
                        # TODO: This also extracts smaller contigs, but eh.
                        for name in names.split(' '):
                            contig_names.add(r.match(name).group(1))
                    for name in contig_names:
                        contig_names_f.write("{}\n".format(name))
                    contig_names_f.flush()

                    contigs_f.close()
                    cmd = f'mfqe --fasta-read-name-lists {contig_names_f.name} --output-fasta-files {contigs_f.name} --output-uncompressed --input-fasta {assembly}'
                    logging.info("Running mfqe to extract {} contigs from {}".format(len(contig_names), assembly))
                    extern.run(cmd)

                # Call prodigal on tempfile
                with tempfile.NamedTemporaryFile(mode='w',suffix='.fasta') as prodigal_f:
                    prodigal_f.close()
                    cmd = f'prodigal -i {contigs_f.name} -a {prodigal_f.name} -c -q -p meta -o /dev/null'
                    logging.info("Running prodigal ..")
                    extern.run(cmd)
            
                    # Diamond against target genome
                    cmd = f'diamond blastp -q {prodigal_f.name} -d {args.target_genome_diamond_db} -p {args.threads}'
                    logging.info("Running diamond against target genome ..")
                    target_diamond_stdout = extern.run(cmd)

                    # Diamond against GTDB
                    cmd = f'diamond blastp -q {prodigal_f.name} -d {args.gtdb_diamond_db} -p {args.threads}'
                    logging.info("Running diamond against GTDB..")
                    gtdb_diamond_stdout = extern.run(cmd)

                    # Parse output
                    # Take the best hit
                    
                    target_df = pd.read_csv(StringIO(target_diamond_stdout),sep='\t',header=None)
                    target_df['db'] = 'target'

                    gtdb_df = pd.read_csv(StringIO(gtdb_diamond_stdout),sep='\t',header=None)
                    gtdb_df['db'] = 'gtdb'

                    df = pd.concat([target_df,gtdb_df])
                    df.columns = ['orf','subject','percent_id','2','3','4','5','6','7','8','evalue','bitscore','db']

                    num_target = 0
                    num_gtdb = 0
                    orf_r = re.compile('(.*)_\d+$')
                    for (orf, d) in df.groupby('orf'):
                        d = d.sort_values('bitscore',ascending=False)
                        d = d.iloc[0]
                        contig = orf_r.match(orf).group(1)
                        if d['db'] == 'target':
                            taxonomy = '-'
                            num_target += 1
                        else:
                            gtdb_accession = d['subject'].split('~')[0]
                            taxonomy = gtdb_metadata[gtdb_metadata['accession'] == gtdb_accession]['gtdb_taxonomy'].values[0]
                            num_gtdb += 1
                        print("\t".join((assembly, contig, d['orf'], d['db'], d['subject'], str(d['percent_id']), str(d['evalue']), str(d['bitscore']), taxonomy)), file=best_hits_output_f)
                    logging.info("Found {} best matching target, and {} best matching gtdb".format(num_target, num_gtdb))

def extract_target_contigs(args):
    # Read in best hits, just going with those that are best matching target
    logging.info("Reading classify contigs output")
    best_hits_df = pd.read_csv(args.classify_contigs_output, sep='\t')
    best_hits_df.columns.values[0] = 'assembly_file'
    best_hits_df.columns.values[1] = 'contig'
    best_hits_df.columns.values[3] = 'hit_type'
    for (assembly, contig_hits) in best_hits_df.groupby('assembly_file'):
        logging.info("Processing {} ..".format(assembly))

        contig_names = set(contig_hits[contig_hits['hit_type']=='target']['contig'].values)
        if len(contig_names) == 0:
            continue
        logging.info("Extracting {} contigs from {}".format(len(contig_names), assembly))

        # Extract contigs into STDOUT
        with tempfile.NamedTemporaryFile(mode='w',suffix='.names') as contigs_f:
            contigs_f.write("\n".join(list(contig_names)))
            contigs_f.flush()

            cmd = f'mfqe --fasta-read-name-lists {contigs_f.name} --output-fasta-files /dev/stdout --output-uncompressed --input-fasta {assembly}'
            contigs = extern.run(cmd)
            # print(cmd)
            # import IPython; IPython.embed()
            # print(contigs)

            for c in StringIO(contigs):
                if c[0]=='>':
                    c = f'>{assembly}:{c[1:]}'
                print(c.strip())
    


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser = argparse.ArgumentParser(parents=[parent_parser])
    subparsers = parser.add_subparsers(title="Sub-commands", dest='subparser_name')

    find_contigs_parser = subparsers.add_parser('find_contigs', help='find contigs containing target OTU sequences')
    find_contigs_parser.add_argument('--query-against-reads',required=True, help="output of 'singlem query' against reads that were used to generate the assembly. Can have used --sequence-type protein or not.")
    find_contigs_parser.add_argument('--sample-list',help="Only consider reads from these samples [default: Consider all]")
    find_contigs_parser.add_argument('--assembly-files',nargs='+',required=True,help="Assembly contig files to search")
    find_contigs_parser.add_argument('--metapackage',required=True)
    find_contigs_parser.add_argument('--threads',default=1,type=int)

    classify_contigs_parser = subparsers.add_parser('classify_contigs', help='classify contigs by protein comparison with GTDB')
    classify_contigs_parser.add_argument('--find-contigs-output',required=True)
    classify_contigs_parser.add_argument('--min-contig-length',default=4000,type=int)
    classify_contigs_parser.add_argument('--gtdb-diamond-db',required=True)
    classify_contigs_parser.add_argument('--target-genome-diamond-db',required=True)
    classify_contigs_parser.add_argument('--threads',default=1,type=int)
    classify_contigs_parser.add_argument('--output-best-hits',required=True)
    classify_contigs_parser.add_argument('--gtdb-bacteria-metadata',required=True)
    classify_contigs_parser.add_argument('--gtdb-archaea-metadata',required=True)

    extract_target_contigs_parser = subparsers.add_parser('extract_target_contigs', help='extract target contigs from assembly')
    extract_target_contigs_parser.add_argument('--classify-contigs-output',required=True)

    args = parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    if args.subparser_name == 'find_contigs':
        find_contigs(args)
    elif args.subparser_name == 'classify_contigs':
        classify_contigs(args)
    elif args.subparser_name == 'extract_target_contigs':
        extract_target_contigs(args)
    else:
        raise Exception("Unknown subcommand")


