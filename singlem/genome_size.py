import logging
import polars as pl

class GenomeSizes:
    @staticmethod
    def calculate_rank_genome_sizes(gtdb_taxonomy_strings, genome_sizes):
        gc = pl.DataFrame({
            "gtdb_taxonomy": gtdb_taxonomy_strings,
            "genome_size": genome_sizes
        })

        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(0).str.strip_chars().alias("kingdom"))
        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(1).str.strip_chars().alias("phylum"))
        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(2).str.strip_chars().alias("class"))
        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(3).str.strip_chars().alias("order"))
        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(4).str.strip_chars().alias("family"))
        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(5).str.strip_chars().alias("genus"))
        gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(6).str.strip_chars().alias("species"))

        # Family level and above are calculated by averaging the genus-wise sizes, to account for oversampling of some genera
        gc = gc.with_columns(pl.col('genome_size').mean().over('genus').alias('genus_wise_genome_size'))

        levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        final = None

        for level in levels:
            if level == 'species':
                taxon_means = gc.select([level, 'genome_size']).group_by(level).mean()
            else:
                taxon_means = gc.group_by('genus').first().select([level, 'genus_wise_genome_size']).group_by(level).mean()
            taxon_means.columns = ['rank', 'genome_size']
            # taxon_mins = gc.select([level, 'genome_size']).group_by(level).min()
            # taxon_mins.columns = ['rank', 'min_genome_size']
            # taxon_maxs = gc.select([level, 'genome_size']).group_by(level).max()
            # taxon_maxs.columns = ['rank', 'max_genome_size']
            # taxon_99pct = gc.select([level, 'genome_size']).group_by(level).quantile(0.99, interpolation='linear')
            # taxon_99pct.columns = ['rank', '99pct_genome_size']

            new_df = taxon_means#.join(taxon_mins, on='rank', how="inner").join(taxon_maxs, on='rank', how="inner").join(taxon_99pct, on='rank', how="inner")

            if final is None:
                final = new_df
            else:
                final = pl.concat([
                    final,
                    new_df
                ])
        logging.info("Calculated genome size statistics for {} taxa".format(len(final)))

        return final

    @staticmethod
    def corrected_genome_size(genome_size, completeness, contamination):
        return float(genome_size) / (1+contamination) / completeness