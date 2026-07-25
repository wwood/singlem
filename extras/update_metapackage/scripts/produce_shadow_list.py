##############################
### produce_shadow_list.py ###
##############################
# Author: Samuel Aroney
# Produce list of non-representative genomes.

import polars as pl

NUM_PER_CLUSTER = 20

bac_metadata = pl.read_csv(snakemake.params.gtdb_bac_metadata, separator="\t", infer_schema_length=10000000)
arc_metadata = pl.read_csv(snakemake.params.gtdb_arc_metadata, separator="\t", infer_schema_length=10000000)
metadata = pl.concat([bac_metadata, arc_metadata])
metadata = metadata.with_columns(
    (pl.col("checkm_completeness") - 5 * pl.col("checkm_contamination")).alias("quality")
).select(["accession", "quality"])

clusters = pl.read_csv(snakemake.params.gtdb_sp_clusters, separator="\t")
clusters = clusters.select(["Representative genome", "Clustered genomes"])
clusters = clusters.with_columns(
    pl.col("Clustered genomes").str.split(",")
).explode("Clustered genomes")
clusters = clusters.filter(pl.col("Representative genome") != pl.col("Clustered genomes"))

shadows = (clusters
    .join(metadata, left_on="Representative genome", right_on="accession")
    .sort("quality", descending=True)
    .group_by("Representative genome", maintain_order=True)
    .head(NUM_PER_CLUSTER)
)

shadows = shadows.with_columns(
    pl.col("Clustered genomes").str.replace_all(r"RS_|GB_", "")
)

# Group into chunks of 1000 and join with commas
shadows = shadows.with_row_index("idx")
shadows = shadows.with_columns((pl.col("idx") // 1000).alias("group"))
shadows = (shadows
    .group_by("group", maintain_order=True)
    .agg(pl.col("Clustered genomes").str.join(","))
)

shadows.select(["group", "Clustered genomes"]).write_csv(
    snakemake.output.genomes, separator="\t", include_header=False
)
