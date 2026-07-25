"""
Test the lyrebird metapackage on three public virome datasets.

Run with:
    snakemake -s 4_test_mpkg.smk --cores 8
"""

import os

output_dir = "test_mpkg_output"
logs_dir = output_dir + "/logs"
PIXI_RUN = "pixi run --manifest-path ~/git/singlem/pixi.toml -e dev"
# Having problems with IO, so use scratch
METAPACKAGE_ARG = "--metapackage /home/woodcrob/s/lyrebird_v0.3.1_phrog_v4.1_metapackage_20250720.smpkg.zb/"

SAMPLES = {
    "ocean_virome":     "ERR2750826",
    "human_gut_virome": "SRR18842432",
    "soil_virome":      "SRR8487022",
}

# All profile outputs that should get krona plots
PROFILE_OUTPUTS = [
    output_dir + "/pipe/ocean_virome.profile.tsv",
    output_dir + "/pipe/human_gut_virome.profile.tsv",
    output_dir + "/pipe/soil_virome.profile.tsv",
    output_dir + "/renew/soil_virome.profile.tsv",
]

# mkdir -p for all output directories
for subdir in ["reads", "pipe", "archive", "renew", "logs/download", "logs/pipe", "logs/pipe_archive", "logs/renew", "logs/krona", "benchmarks/download", "benchmarks/pipe", "benchmarks/pipe_archive", "benchmarks/renew", "benchmarks/krona"]:
    os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)

rule all:
    input:
        # Standard pipe -p for ocean and human gut
        expand(output_dir + "/pipe/{sample}.profile.tsv", sample=["ocean_virome", "human_gut_virome"]),
        # Soil: pipe -p (direct)
        output_dir + "/pipe/soil_virome.profile.tsv",
        # Soil: pipe --archive-otu-table then renew -p
        output_dir + "/renew/soil_virome.profile.tsv",
        # Krona plots for all profiles
        [p + ".krona.html" for p in PROFILE_OUTPUTS],

rule download:
    output:
        r1 = output_dir + "/reads/{sample}_1.fastq.gz",
        r2 = output_dir + "/reads/{sample}_2.fastq.gz",
    params:
        accession = lambda wc: SAMPLES[wc.sample],
        outdir = output_dir + "/reads",
    log:
        output_dir + "/logs/download/{sample}.log",
    benchmark:
        output_dir + "/benchmarks/download/{sample}.txt"
    resources:
        mem_mb = 4 * 1024,
        runtime = 4 * 60,
    shell:
        # Current released version of kingfisher fails for prefetch due to deprecated --output-file, so use bleeding edge.
        """
        (pixi run --manifest-path ~/git/kingfisher/pixi.toml kingfisher get \
            -r {params.accession} \
            -m ena-ftp prefetch \
            --output-format-possibilities fastq.gz \
            --output-directory {params.outdir} && \
        mv {params.outdir}/{params.accession}_1.fastq.gz {output.r1} && \
        mv {params.outdir}/{params.accession}_2.fastq.gz {output.r2}) &> {log}
        """

rule pipe_profile:
    input:
        r1 = output_dir + "/reads/{sample}_1.fastq.gz",
        r2 = output_dir + "/reads/{sample}_2.fastq.gz",
    output:
        profile = output_dir + "/pipe/{sample}.profile.tsv",
        archive = output_dir + "/pipe/{sample}.archive.otu_table.tsv",
    log:
        output_dir + "/logs/pipe/{sample}.log",
    benchmark:
        output_dir + "/benchmarks/pipe/{sample}.txt"
    threads: 8
    resources:
        mem_mb = 16 * 1024,
        runtime = 12 * 60,
    shell:
        """
        ({PIXI_RUN} lyrebird pipe \
            {METAPACKAGE_ARG} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {output.profile} \
            --archive-otu-table {output.archive} \
            --threads {threads}) &> {log}
        """

rule pipe_archive:
    """Run lyrebird pipe --archive-otu-table on soil virome."""
    input:
        r1 = output_dir + "/reads/soil_virome_1.fastq.gz",
        r2 = output_dir + "/reads/soil_virome_2.fastq.gz",
    output:
        archive = output_dir + "/archive/soil_virome.archive.otu_table.tsv",
    log:
        output_dir + "/logs/pipe_archive/soil_virome.log",
    benchmark:
        output_dir + "/benchmarks/pipe_archive/soil_virome.txt"
    threads: 8
    resources:
        mem_mb = 16 * 1024,
        runtime = 12 * 60,
    shell:
        """
        ({PIXI_RUN} lyrebird pipe \
            {METAPACKAGE_ARG} \
            -1 {input.r1} \
            -2 {input.r2} \
            --archive-otu-table {output.archive} \
            --no-assign-taxonomy \
            --threads {threads}) &> {log}
        """

rule krona:
    """Generate a Krona plot from a taxonomic profile."""
    input:
        profile = output_dir + "/{path}.profile.tsv",
    output:
        krona = output_dir + "/{path}.profile.tsv.krona.html",
    log:
        output_dir + "/logs/krona/{path}.log",
    benchmark:
        output_dir + "/benchmarks/krona/{path}.txt"
    resources:
        mem_mb = 4 * 1024,
        runtime = 30,
    shell:
        """
        ({PIXI_RUN} singlem summarise \
            --input-taxonomic-profiles {input.profile} \
            --output-taxonomic-profile-krona {output.krona}) &> {log}
        """

rule renew_profile:
    """Run lyrebird renew -p on the archive OTU table from the soil virome."""
    input:
        archive = output_dir + "/archive/soil_virome.archive.otu_table.tsv",
    output:
        profile = output_dir + "/renew/soil_virome.profile.tsv",
    log:
        output_dir + "/logs/renew/soil_virome.log",
    benchmark:
        output_dir + "/benchmarks/renew/soil_virome.txt"
    threads: 8
    resources:
        mem_mb = 16 * 1024,
        runtime = 12 * 60,
    shell:
        """
        ({PIXI_RUN} lyrebird renew \
            {METAPACKAGE_ARG} \
            --input-archive-otu-table {input.archive} \
            -p {output.profile} \
            --threads {threads}) &> {log}
        """
