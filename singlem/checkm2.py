import polars as pl
import logging
from dataclasses import dataclass

@dataclass
class CheckM2Stats:
    completeness: float
    contamination: float

@dataclass
class CheckM2MoreStats:
    completeness: float
    contamination: float
    genome_size: int

class CheckM2:
    def __init__(self, quality_file):
        self.quality_file = quality_file
        self.qualities = pl.read_csv(self.quality_file, separator='\t')
        logging.info("Read in {} genome qualities".format(self.qualities.shape[0]))

    def genomes_of_sufficient_quality(self, min_completeness, max_contamination):
        return list(
            self.qualities.filter((pl.col("Completeness") >= min_completeness) &
                                  (pl.col("Contamination") <= max_contamination)).select('Name').get_columns()[0])

    # Implement "in"
    def __contains__(self, item):
        return self.qualities.filter(pl.col("Name") == item).shape[0] > 0

    def names(self):
        '''Return a list of all genome names read in'''
        return list(self.qualities['Name'].to_list())

    def get_stats(self, genome_name):
        found = self.qualities.filter(pl.col("Name") == genome_name)
        if found.shape[0] != 1:
            raise Exception("Expected to find 1 genome, found {}".format(found.shape[0]))
        
        for row in found.rows(named=True):
            return CheckM2Stats(
                completeness=row['Completeness'] / 100.,
                contamination=row['Contamination'] / 100.
            )

    def get_all_stats(self):
        stats = {}
        for row in self.qualities.rows(named=True):
            stats[row['Name']] = CheckM2MoreStats(
                completeness=row['Completeness'] / 100.,
                contamination=row['Contamination'] / 100.,
                genome_size=row['Genome_Size'],
            )
        return stats