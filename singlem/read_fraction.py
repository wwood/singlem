import pandas as pd
import polars as pl
import logging
import json
import os
import sys
from dataclasses import dataclass
from heapq import heappush, heappushpop

import extern

from .condense import CondensedCommunityProfile
from .metapackage import Metapackage
from .utils import FastaNameToSampleName

@dataclass
class GenomeSizeStruct:
    mean: float

class ReadFractionEstimator:
    def calculate_and_report_read_fraction(self, **kwargs):
        # input_profile = args.input_profile,
        # metagenome_sizes = args.input_metagenome_sizes,
        # taxonomic_genome_lengths_file = args.taxonomic_genome_lengths_file
        input_profile = kwargs.pop('input_profile')
        metagenome_sizes = kwargs.pop('metagenome_sizes')
        forward_read_files = kwargs.pop('forward_read_files')
        reverse_read_files = kwargs.pop('reverse_read_files')
        taxonomic_genome_lengths_file = kwargs.pop('taxonomic_genome_lengths_file')
        metapackage = kwargs.pop('metapackage')
        accept_missing_samples = kwargs.pop('accept_missing_samples')
        output_tsv = kwargs.pop('output_tsv')
        output_per_taxon_read_fractions = kwargs.pop('output_per_taxon_read_fractions')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        output_fh = open(output_tsv, 'w') if output_tsv else sys.stdout

        if output_per_taxon_read_fractions:
            output_per_taxon_read_fractions_fh = open(output_per_taxon_read_fractions, 'w')

        # Grab the genome length data
        if taxonomic_genome_lengths_file:
            taxonomic_genome_lengths_df = pd.read_csv(taxonomic_genome_lengths_file, sep='\t')
        else:
            if metapackage:
                mpkg = Metapackage.acquire(metapackage)
            else:
                mpkg = Metapackage.acquire_default()
            if mpkg.version < 4:
                raise Exception("Taxonomic genome lengths files are only included in v4+ metapackages.")
            taxonomic_genome_lengths_df = mpkg.taxon_genome_lengths()

        # Read in the taxonomic genome lengths
        if 'rank' not in taxonomic_genome_lengths_df.columns:
            raise Exception("Taxonomic genome lengths file must have a 'rank' column.")
        if 'genome_size' not in taxonomic_genome_lengths_df.columns:
            raise Exception("Taxonomic genome lengths file must have a 'genome_size' column.")
        taxonomic_genome_lengths = {}
        for _, row in taxonomic_genome_lengths_df.iterrows():
            taxonomic_genome_lengths[row['rank']] = GenomeSizeStruct(row['genome_size'])
        logging.info("Read taxonomic genome lengths for %i rank(s)." % len(taxonomic_genome_lengths))

        # Read in the metagenome sizes
        if metagenome_sizes:
            if forward_read_files or reverse_read_files:
                raise Exception("Cannot specify both a metagenome sizes file and read files.")
            metagenome_sizes_df = pd.read_csv(metagenome_sizes, sep='\t')
            if 'sample' not in metagenome_sizes_df.columns:
                raise Exception("Metagenome sizes file must have a 'sample' column.")
            if 'num_bases' not in metagenome_sizes_df.columns:
                raise Exception("Metagenome sizes file must have a 'num_bases' column.")
            metagenome_sizes = {}
            for _, row in metagenome_sizes_df.iterrows():
                metagenome_sizes[row['sample']] = float(row['num_bases'])
            logging.info("Read metagenome sizes for %i sample(s)" % len(metagenome_sizes))
        else:
            if not forward_read_files:
                raise Exception("Must specify either a metagenome sizes file or read files.")
            metagenome_sizes = self._get_stems_and_read_files(forward_read_files, reverse_read_files)
            
        print("\t".join([
            "sample",
            "bacterial_archaeal_bases",
            "metagenome_size",
            "read_fraction",
            "average_bacterial_archaeal_genome_size",
            "warning",
            "domain_relative_abundance",
            "phylum_relative_abundance",
            "class_relative_abundance",
            "order_relative_abundance",
            "family_relative_abundance",
            "genus_relative_abundance",
            "species_relative_abundance"]),
            file=output_fh)

        if output_per_taxon_read_fractions:
            print("sample\ttaxonomy\tbase_contribution", file=output_per_taxon_read_fractions_fh)
            
        num_samples = 0
        # Iterate through the input profile, calculating the read fraction for each sample
        with open(input_profile) as f:
            for profile in CondensedCommunityProfile.each_sample_wise(f):
                sample = profile.sample
                if sample not in metagenome_sizes:
                    if accept_missing_samples:
                        logging.warning("Sample '%s' in profile not found in metagenome sizes file. Skipping." % sample)
                        continue
                    else:
                        raise Exception("Sample '%s' in profile not found in metagenome sizes file." % sample)
                metagenome_size = metagenome_sizes[sample]

                # Calculate the read fraction. This is the coverage of each taxon multiplied by its average genome size.
                account = 0
                # Priority queue with highest abundances first
                num_unknown_taxa_accounted_for = 3
                @dataclass(order=True)
                class PrioritizedTaxon:
                    priority: float  # Negative of relabund to make it a max queue
                    taxon: str

                highest_unknown_taxon_names = []
                min_unknown_taxon_priority = 0  # To bound the size of the heapq

                average_genome_size_numerator = 0
                average_genome_size_denominator = 0

                for node in profile.breadth_first_iter():
                    taxonomy = node.word
                    if taxonomy == 'Root':
                        continue  # More likely false positive hits, I guess.
                    if taxonomy not in taxonomic_genome_lengths:
                        raise Exception("Taxonomy '%s' in profile not found in taxonomic genome lengths file." % taxonomy)
                    genome_size = taxonomic_genome_lengths[taxonomy].mean
                    contribution = node.coverage * genome_size
                    account += contribution

                    average_genome_size_numerator += node.coverage * genome_size
                    average_genome_size_denominator += node.coverage

                    if '__' not in taxonomy or node.calculate_level() > 7:
                        raise Exception("It appears that a condensed profile with a non-standard taxonomy was used. This is not supported for prokaryotic_fraction (microbial_fraction), since prokaryotic_fraction relies on knowing whether the taxon has been assigned to the species level or not.")
                    if 's__' not in taxonomy:
                        if node.coverage > min_unknown_taxon_priority:
                            if len(highest_unknown_taxon_names) < num_unknown_taxa_accounted_for:
                                # At start when queue is not full, don't remove anything
                                heappush(highest_unknown_taxon_names, PrioritizedTaxon(-node.coverage, taxonomy))
                            else:
                                heappushpop(highest_unknown_taxon_names, PrioritizedTaxon(-node.coverage, taxonomy))
                                # minus minus since we want to remove the lowest priority
                                min_unknown_taxon_priority = -highest_unknown_taxon_names[-1].priority

                    if contribution > 0 and output_per_taxon_read_fractions:
                        print("%s\t%s\t%s" % (sample, taxonomy, contribution),
                            file=output_per_taxon_read_fractions_fh)
                
                highest_unknown_taxa = highest_unknown_taxon_names
                highest_unknown_taxa_sum = sum([-x.priority * taxonomic_genome_lengths[x.taxon].mean for x in highest_unknown_taxa])
                doubled_account = account + highest_unknown_taxa_sum
                halved_account = account - (highest_unknown_taxa_sum / 2)
                warning_threshold = 0.02
                if (doubled_account / metagenome_size - account / metagenome_size) > warning_threshold or \
                    (account / metagenome_size - halved_account / metagenome_size) > warning_threshold:
                    warning = "WARNING: The most abundant taxons not assigned to the species level account for a large fraction of the total estimated read fraction. This may mean that the read_fraction estimate is inaccurate."
                else:
                    warning = ""

                final_estimate = account / metagenome_size * 100
                if final_estimate > 100:
                    final_estimate = 100

                # Average genome size is the weighted average of the genome sizes of the taxa in the profile
                if average_genome_size_denominator == 0:
                    average_genome_size = 0
                else:
                    average_genome_size = average_genome_size_numerator / average_genome_size_denominator

                per_level_coverage = profile.taxonomic_level_coverage_table(assume_8_levels=True)

                print("%s\t%s\t%s\t%0.2f\t%s\t%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f" % (
                    sample,
                    round(account),
                    metagenome_size,
                    final_estimate,
                    round(average_genome_size),
                    warning,
                    per_level_coverage.filter(pl.col('level') == 'domain')['relative_abundance'][0],
                    per_level_coverage.filter(pl.col('level') == 'phylum')['relative_abundance'][0],
                    per_level_coverage.filter(pl.col('level') == 'class')['relative_abundance'][0],
                    per_level_coverage.filter(pl.col('level') == 'order')['relative_abundance'][0],
                    per_level_coverage.filter(pl.col('level') == 'family')['relative_abundance'][0],
                    per_level_coverage.filter(pl.col('level') == 'genus')['relative_abundance'][0],
                    per_level_coverage.filter(pl.col('level') == 'species')['relative_abundance'][0]),
                    file=output_fh)

                num_samples += 1
        logging.info("Calculated read fractions for %d samples." % num_samples)

        if output_tsv: output_fh.close()
        if output_per_taxon_read_fractions: output_per_taxon_read_fractions_fh.close()

        logging.info("Finished.")

    def _get_stems_and_read_files(self, forward_read_files, reverse_read_files):
        stems_to_read_files = {}
        for i in range(len(forward_read_files)):
            if not os.path.exists(forward_read_files[i]):
                raise Exception("Forward read file '%s' does not exist." % forward_read_files[i])
            if reverse_read_files:
                if not os.path.exists(reverse_read_files[i]):
                    raise Exception("Reverse read file '%s' does not exist." % reverse_read_files[i])
            f = forward_read_files[i]
            stem = FastaNameToSampleName.fasta_to_name(f)

            to_add = [f]
            if reverse_read_files:
                to_add.append(reverse_read_files[i])
            if stem in stems_to_read_files:
                raise Exception("Found duplicate stem '%s' in forward read files." % stem)
            stems_to_read_files[stem] = to_add

        logging.debug("Found stems and read files: %s" % stems_to_read_files)
        logging.info("Analysing %i metagenome(s) .." % len(stems_to_read_files))
        return SmafaCountedMetagenomeSizes(stems_to_read_files)


class SmafaCountedMetagenomeSizes:
    def __init__(self, stems_to_read_files):
        self.stems_to_read_files = stems_to_read_files

    def __contains__(self, stem):
        # Somewhat of a hack but allow fuzziness for _1
        return stem in self.stems_to_read_files or stem + '_1' in self.stems_to_read_files

    def __getitem__(self, stem):
        if stem not in self:
            raise Exception("Stem '%s' not found in input metagenome set." % stem)
        total_base_count = 0
        logging.info("Counting bases in sample '%s' .." % stem)
        try:
            read_files = self.stems_to_read_files[stem]
        except KeyError:
            read_files = self.stems_to_read_files[stem + '_1']
        j = extern.run('smafa count -i %s' % ' '.join(read_files))
        logging.debug("Found JSON response from smafa count: %s" % j)
        j2 = json.loads(j)
        
        # [{"path":"/dev/fd/63","num_reads":1,"num_bases":3},{"path":"/dev/fd/62","num_reads":1,"num_bases":2}]
        for read_file in j2:
            total_base_count += int(read_file['num_bases'])
        logging.info("Total base count for sample '%s' is %.2f Gbp" % (stem, total_base_count / 1_000_000_000))

        return total_base_count