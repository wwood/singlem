#!/usr/bin/env python

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2017"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

import argparse
import logging
import os
import sys
import csv
import tempfile
import extern
import re

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

import singlem.pipe as pipe
from singlem.querier import Querier, QueryInputSequence
from singlem.sequence_database import SequenceDatabase
from singlem.sequence_classes import SeqReader
from singlem.archive_otu_table import ArchiveOtuTable
from singlem.checkm import CheckM


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=singlem.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--contigs', required=True)
    parser.add_argument('--checkm_csv', required=True)
    parser.add_argument('--bin_files', nargs='+', required=True)
    parser.add_argument('--db', required=True)

    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--collapse_sra_coupled', action='store_true', help='Remove _1 and _2 from sample names')
    parser.add_argument('--bin_file_extension', default="fna")
    parser.add_argument('--samples_to_ignore', nargs='+', default=[])
    parser.add_argument('--samples_to_pick', type=int, default=5)
    parser.add_argument('--samples_per_otu', type=int, default=1)

    parser.add_argument('--singlem_on_contigs_archive_otu_table')

    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    samples_to_ignore = set(args.samples_to_ignore)
    num_samples_to_pick = args.samples_to_pick

    # Read checkm CSV file
    with open(args.checkm_csv) as f:
        quality_bins = CheckM.read_checkm_stats(f, 80, 5)
    logging.info("Found %i bin(s) of sufficient quality according to CheckM statistics" % len(quality_bins))

    # Read in contig lengths - no interested in contigs < 2.5kb a la metabat
    contig_to_length = {}
    with open(args.contigs) as io:
        for contig_name, seq, _ in SeqReader().readfq(io):
            if contig_name in contig_to_length:
                raise Exception("Duplicate contig name in FASTA file detected: %s" % contig_name)
            contig_to_length[contig_name] = len(seq)
    logging.info("Read in contig lengths for %i contigs" % len(contig_to_length))

    # Read bin files, recording which contigs made it into each bin
    bin_names = set()
    contig_to_bin = {}
    bin_extension_regex = re.compile(r'\.%s$' % args.bin_file_extension)
    for bin_file in args.bin_files:
        base = os.path.basename(bin_file)
        if not bin_extension_regex.search(base):
            raise Exception("Bin name %s did not match regex" % bin_file)
        bin_name = bin_extension_regex.sub('',base)
        logging.debug("Using bin name %s" % bin_name)
        if bin_name in bin_names:
            raise Exception("Duplicate bin names detected: %s" % bin_name)
        if bin_name in quality_bins:
            logging.debug("Bin %s is of sufficient quality, so reading in contig names" % bin_name)
            with open(bin_file) as f:
                for contig_name, _, _ in SeqReader().readfq(f):
                    contig_to_bin[contig_name] = bin_name
        else:
            logging.debug("Bin %s not of sufficient quality, so ignoring contained contigs" % bin_name)
    logging.info("Read in binning information for %i contigs" % len(contig_to_bin))

    # Run SingleM on contigs

    if args.singlem_on_contigs_archive_otu_table:
        with open(args.singlem_on_contigs_archive_otu_table) as a:
            otus = list(ArchiveOtuTable.read(a))
    else:
        with tempfile.NamedTemporaryFile(prefix="singlem-binning") as tf:
            # TODO: Run this internally rather than through extern. Need to have better default
            # options in Pipe.run
            logging.info("Running singlem pipe on the contigs..")
            json_file = 'singlem_on_contigs.json'#tf.name,
            extern.run("singlem pipe --sequences '%s' --archive_otu_table '%s' --threads %i "
                       "--output_extras --no_assign_taxonomy" %(
                           args.contigs, json_file, args.threads))
            with open(json_file) as a:
                otus = list(ArchiveOtuTable.read(a))
    logging.info("Finished running singlem pipe, found %i OTUs amongst the different genes" % len(otus))

    # If there is any more genomes
    # Which 2.5kb contigs have OTUs that are not binned?
    unbinned_contigs = set()
    unbinned_otus = []
    num_otus_binned = 0
    num_unbinned_otus_in_short_contigs = 0
    num_unbinned_otus_in_long_contigs = 0
    length_cutoff = 2500
    for otu in otus:
        contig_names_field = otu.fields.index('read_names') # There should be a better way to do this.
        any_unbinned = False
        for contig_name in otu.data[contig_names_field]:
            if contig_name not in contig_to_length:
                raise Exception("Unexpected contig found in SingleM file: %s" % contig_name)
            if contig_name in contig_to_bin:
                num_otus_binned += 1
            elif contig_to_length[contig_name] >= length_cutoff:
                num_unbinned_otus_in_long_contigs += 1
                unbinned_contigs.add(contig_name)
                if any_unbinned == False:
                    unbinned_otus.append(otu)
                    any_unbinned = True
            elif contig_to_length[contig_name] < length_cutoff:
                num_unbinned_otus_in_short_contigs += 1
            else: raise Exception("Programming error")
    logging.info("Found %i OTUs in good bins, %i insufficiently binned OTUs in "
                 "contigs >= %ibp, and %i insufficiently binned OTUs in shorter contigs" % (
                     num_otus_binned, num_unbinned_otus_in_long_contigs,
                     length_cutoff, num_unbinned_otus_in_short_contigs))
    if len(unbinned_otus) == 0:
        logging.info("All OTUs detected in sufficiently long contigs are already in "
                     "quality bins. Congratulations. SingleM cannot really help you any further.")
    else:
        query_results = []
        unbinned_seqs = set()
        for otu in unbinned_otus:
            unbinned_seqs.add(otu.sequence)
        logging.info("Before picking samples, found %i OTU sequences" % len(unbinned_seqs))
        picked_samples = set()
        num_samples_per_otu = {}
        for seq in unbinned_seqs:
            num_samples_per_otu[seq] = 0

        while len(picked_samples) < num_samples_to_pick:
            sample_to_seqs_detected = {}
            queries = [QueryInputSequence(sequence,sequence) for sequence in unbinned_seqs]

            db = SequenceDatabase.acquire(args.db)
            max_divergence = 0
            max_target_seqs = 1000 # Not actually used since SQLite is used, not BLAST
            num_threads = 1
            query_result = Querier().query_with_queries(
                queries, db, max_target_seqs, max_divergence, num_threads)

            for result in query_result:
                sample = result.subject.sample_name
                if args.collapse_sra_coupled:
                    sample_name = sample.replace('_1','').replace('_2','')
                else:
                    sample_name = sample
                if sample_name not in samples_to_ignore:
                    seq = result.query.sequence
                    try:
                        sample_to_seqs_detected[sample_name].add(seq)
                    except:
                        sample_to_seqs_detected[sample_name] = set([seq])

            max_sample = None
            max_len = 0
            for sample, seqs in sample_to_seqs_detected.items():
                if sample not in picked_samples and (
                        len(seqs) > max_len or (
                            # Choose lowest lexical for reproducibility.
                            len(seqs) == max_len and sample < max_sample)):
                    logging.debug("%s has %i unbinned seqs" % (sample, len(seqs)))
                    max_sample = sample
                    max_len = len(seqs)
            if max_sample:
                picked_samples.add(max_sample)
                print max_sample
                for seq in sample_to_seqs_detected[max_sample]:
                    num_samples_per_otu[seq] += 1
                    if num_samples_per_otu[seq] >= args.samples_per_otu:
                        unbinned_seqs.remove(seq)
            else:
                logging.info("Found no more samples containing unbinned reads not already accounted for")
                break


