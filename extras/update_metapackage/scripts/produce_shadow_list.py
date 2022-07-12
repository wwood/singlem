##############################
### produce_shadow_list.py ###
##############################
# Author: Samuel Aroney
# Produce list of non-representative genomes.

import pandas as pd

NUM_PER_CLUSTER = 20

bac_metadata = pd.read_csv(snakemake.params.gtdb_bac_metadata, sep = "\t")
arc_metadata = pd.read_csv(snakemake.params.gtdb_arc_metadata, sep = "\t")
metadata = pd.concat([bac_metadata, arc_metadata])
metadata["quality"] = metadata["checkm_completeness"] - 5 * metadata["checkm_contamination"]
metadata = metadata[["accession", "quality"]].set_index("accession")

clusters = pd.read_csv(snakemake.params.gtdb_sp_clusters, sep = "\t")
clusters = clusters[["Representative genome", "Clustered genomes"]]
clusters["Clustered genomes"] = clusters["Clustered genomes"].str.split(",")
clusters = clusters.explode("Clustered genomes")
clusters = clusters[clusters["Representative genome"] != clusters["Clustered genomes"]]

shadows = (clusters
    .join(metadata, on = "Representative genome")
    .sort_values("quality", ascending=False)
    .groupby("Representative genome")
    .head(NUM_PER_CLUSTER)
    )

shadows[["Clustered genomes"]].to_csv(snakemake.output.genomes, sep="\t", index=False, header=False)
