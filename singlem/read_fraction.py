import pandas as pd
import logging

from .condense import CondensedCommunityProfile
from .metapackage import Metapackage

class ReadFractionEstimator:
    def calculate_and_report_read_fraction(self, **kwargs):
        # input_profile = args.input_profile,
        # metagenome_sizes = args.input_metagenome_sizes,
        # taxonomic_genome_lengths_file = args.taxonomic_genome_lengths_file
        input_profile = kwargs.pop('input_profile')
        metagenome_sizes = kwargs.pop('metagenome_sizes')
        taxonomic_genome_lengths_file = kwargs.pop('taxonomic_genome_lengths_file')
        metapackage = kwargs.pop('metapackage')
        accept_missing_samples = kwargs.pop('accept_missing_samples')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

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
        for index, row in taxonomic_genome_lengths_df.iterrows():
            taxonomic_genome_lengths[row['rank']] = row['genome_size']
        logging.info("Read taxonomic genome lengths for %i rank(s)." % len(taxonomic_genome_lengths))

        # Read in the metagenome sizes
        metagenome_sizes_df = pd.read_csv(metagenome_sizes, sep='\t')
        if 'sample' not in metagenome_sizes_df.columns:
            raise Exception("Metagenome sizes file must have a 'sample' column.")
        if 'num_bases' not in metagenome_sizes_df.columns:
            raise Exception("Metagenome sizes file must have a 'num_bases' column.")
        metagenome_sizes = {}
        for index, row in metagenome_sizes_df.iterrows():
            metagenome_sizes[row['sample']] = float(row['num_bases'])
        logging.info("Read metagenome sizes for %i sample(s)" % len(metagenome_sizes))

        # Iterate through the input profile, calculating the read fraction for each sample
        read_fractions = {}
        print("sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction")
        num_samples = 0
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
                for node in profile.breadth_first_iter():
                    taxonomy = node.word
                    if taxonomy == 'Root': continue # More likely false positive hits, I guess.
                    if taxonomy not in taxonomic_genome_lengths:
                        raise Exception("Taxonomy '%s' in profile not found in taxonomic genome lengths file." % taxonomy)
                    contribution = node.coverage * taxonomic_genome_lengths[taxonomy]
                    account += contribution
                
                print("%s\t%s\t%s\t%0.2f%%" % (sample, account, metagenome_size, account / metagenome_size * 100))
                num_samples += 1
        logging.info("Calculated read fractions for %d samples." % num_samples)
        logging.info("Finished.")

