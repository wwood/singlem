
import polars as pl
import logging

class CheckM2:
    def __init__(self, quality_file='\t'):
        self.quality_file = quality_file
        self.qualities = pl.read_csv(self.quality_file, separator='\t')
        logging.info("Read in {} genome qualities".format(self.qualities.shape[0]))

    def genomes_of_sufficient_quality(self, min_completeness, max_contamination):
        return list(self.qualities.filter(
            (pl.col("Completeness") >= min_completeness) &
            (pl.col("Contamination") <= max_contamination)
        ).select('Name').get_columns()[0])
