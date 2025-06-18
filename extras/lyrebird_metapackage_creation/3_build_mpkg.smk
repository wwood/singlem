#############
### Setup ###
#############
"""
Steps:
Create and activate base environment (env.yml)
Update config.yaml
Run `snakemake --cores 64 --use-conda --scheduler greedy`
"""

import pandas as pd
import os
# import shutil

output_dir = config["output_dir"]
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
logs_dir = output_dir + "/logs"
if not os.path.exists(logs_dir):
    os.mkdir(logs_dir)
scripts_dir = output_dir + "/qsub_scripts"
hmms_and_names = pd.read_csv(output_dir + "/hmms_and_names_roundrobin.tsv", sep="\t").set_index("gene", drop=False)

rule all:
    input:
        output_dir + "/metapackage.done"

rule chainsaw:
    input:
        spkg = output_dir + "/regenerate/{spkg}.spkg",
    output:
        directory(output_dir + "/chainsaw/{spkg}.spkg"),
    log:
        log = output_dir + "/logs/chainsaw/{spkg}.log"
    resources:
        mem_mb = 4 * 1024,
        runtime = 1 * 60
    conda:
        "envs/singlem.yml"
    shell:
        "singlem chainsaw --input-singlem-package {input.spkg} --output-singlem-package {output} 2> {log}"

rule create_draft_metapackage:
    input:
        packages = expand(output_dir + "/chainsaw/{spkg}.spkg", spkg=hmms_and_names.index),
    output:
        metapackage = directory(output_dir + "/draft_metapackage.smpkg"),
        done = output_dir + "/draft_metapackage.done"
    threads:
        workflow.cores
    log:
        logs_dir + "/draft_metapackage.log"
    resources:
        mem_mb = 128 * 1024,
        runtime = 24 * 2 * 60
    conda:
        "singlem"
    shell:
        "singlem metapackage "
        "--singlem-packages {input.packages} "
        "--no-nucleotide-sdb "
        "--no-taxon-genome-lengths "
        "--makeidx-sensitivity-params '--sensitive ' "
        "--metapackage {output.metapackage} "
        "--threads {threads} "
        "&> {log}; "
        "touch {output.done}"

rule Lyrebird_transcripts:
    input:
        dir = config["viral_genome_fna_reps"],
        metapackage = output_dir + "/draft_metapackage.smpkg",
        done = output_dir + "/draft_metapackage.done"
    output:
        dir = directory(output_dir + "/transcripts"),
        touch = output_dir + "/transcripts.done"
    params:
        logs = logs_dir + "/transcripts",
    conda:
        "envs/singlem.yml"
    threads:
        workflow.cores
    resources:
        mem_mb = 16 * 1024,
        runtime = 12 * 60
    shell:
        "mkdir -p {params.logs} "
        "&& mkdir -p {output.dir} "
        "&& cat <( find {input.dir} -name '*.fna' ) "
        "| parallel --will-cite --eta -j {threads} "
        "singlem pipe "
        "--forward {{}} "
        "--metapackage {input.metapackage} "
        "--otu-table {output.dir}/{{/.}}.otu_table.tsv "
        "--no-assign-taxonomy "
        "'&>' {params.logs}/{{/.}}.log; "
        "touch {output.touch}"

rule assign_taxonomy:
    input:
        touch = output_dir + "/transcripts.done",
        metapackage = output_dir + "/draft_metapackage.smpkg",
    output:
        output_dir + "/assign_taxonomy/transcripts.otu_table.tsv"
    params:
        input_dir = output_dir + "/transcripts",
        viral_taxonomy = config["viral_tax"],
    conda:
        "envs/singlem.yml"
    resources:
        mem_mb = 8 * 1024,
        runtime = 8 * 60
    script:
        "scripts/assign_viral_taxonomy.py"

rule make_sdb:
    input:
        output_dir + "/assign_taxonomy/transcripts.otu_table.tsv"
    output:
        directory(output_dir + "/taxonomy/transcripts.sdb")
    threads: config["max_threads"]
    log:
        logs_dir + "/taxonomy/makedb.log"
    conda:
        "envs/singlem.yml"
    resources:
        mem_mb = 8 * 1024,
        runtime = 2 * 60
    shell:
        "singlem makedb "
        "--otu-table {input} "
        "--db {output} "
        "--threads {threads} "
        "&> {log}"

rule create_Lyrebird_metapackage:
    input:
        packages = expand(output_dir + "/chainsaw/{spkg}.spkg", spkg=hmms_and_names.index),
        sdb = output_dir + "/taxonomy/transcripts.sdb",
    output:
        metapackage = directory(output_dir + "/metapackage/" + config["output_metapackage"]),
        done = output_dir + "/metapackage.done"
    log:
        logs_dir + "/metapackage.log"
    threads:
        workflow.cores
    resources:
        mem_mb = 128 * 1024
        runtime = 24 * 2 * 60
    conda:
        "singlem"
    shell:
        "singlem metapackage "
        "--singlem-packages {input.packages} "
        "--nucleotide-sdb {input.sdb} "
        "--no-taxon-genome-lengths "
        "--makeidx-sensitivity-params '--sensitive ' "
        "--calculate-average-num-genes-per-species "
        "--metapackage {output.metapackage} "
        "--threads {threads} "
        "&> {log}; "
        "touch {output.done}"