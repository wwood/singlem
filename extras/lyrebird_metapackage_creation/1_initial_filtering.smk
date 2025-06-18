#############
### Setup ###
#############
"""
Steps:
Create and activate base environment (env.yml)
Update config.yaml
Run `snakemake --cores 64 --use-conda`
"""

localrules: all, hmmsearch_viral, resolve_conflicts

import pandas as pd
import os

output_dir = config["output_dir"]
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
logs_dir = output_dir + "/logs"
if not os.path.exists(logs_dir):
    os.mkdir(logs_dir)
scripts_dir = output_dir + "/qsub_scripts"
if not os.path.exists(scripts_dir):
    os.mkdir(scripts_dir)

hmm_file = output_dir + "/all_hmms.hmm"
if not os.path.exists(hmm_file):
    with open(hmm_file, "w+") as outfile:
        print("Concatenating hmms into a single file")
        with open(config["hmms_and_names"], "r") as hmm_info:
            for line in hmm_info:
                if not line.startswith("gene"):
                    __, __, hmm_filepath, __ = line.strip().split("\t")
                    with open(hmm_filepath, "r") as infile:
                        outfile.write(infile.read())

hmms_and_names = pd.read_csv(config["hmms_and_names"], sep="\t").set_index("gene", drop=False)

gtdb_proviruses = 'gtdb_proviruses_r214.tsv'
gtdb_refseq_to_assembly = 'assembly_to_refseq.txt'

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
logs_dir = os.path.join(output_dir, "logs")

# if not "max_threads" in config:
#     config["max_threads"] = 64

rule all:
    input:
        output_dir + "/resolve_conflicts.done"

###################################
# Initial GraftM package creation #
###################################
rule hmmsearch_viral:
    input:
        hmm = hmm_file,
        genome_proteins = config["viral_faa_list"]
    output:
        script_dir = temp(directory(scripts_dir + "/hmmsearch_viral")),
        touch = output_dir + "/hmmsearch_viral.done"
    params:
        output_dir = output_dir + "/hmmsearch_viral",
    threads: workflow.cores
    resources:
        mem_mb = 128 * 1024,
        runtime = 24 * 2 * 60
    log:
        logs_dir + "/hmmsearch.log"
    conda:
        "envs/hmmsearch.yml"
    script:
        "scripts/hmmsearch.py"

rule get_matches_viral:
    input:
        touch = output_dir + "/hmmsearch_viral.done"
    output:
        touch = output_dir + "/get_matches_viral.done"
    params:
        proviruses = False,
        hmms_and_names = config["hmms_and_names"],
        hmmsearch_directory = output_dir + "/hmmsearch_viral",
        output_dir = output_dir + "/get_matches_viral",
        logs_dir = logs_dir + "/get_matches_viral"
    threads: workflow.cores
    resources:
        mem_mb = 128 * 1024,
        runtime = 8 * 60
    log:
        logs_dir + "/get_matches_viral.log"
    conda:
        "envs/get_matches.yml"
    script:
        "scripts/get_matches_all.py"

rule resolve_conflicts:
    # take all sequences, if multiple spkgs have the same sequence, take the one with the highest coverage, return list of non-conflict spkgs
    input:
        touch = output_dir + "/get_matches_viral.done"
    output:
        touch = output_dir + "/resolve_conflicts.done",
        hmms_and_names_noconflict = output_dir + "/hmms_and_names_noconflict.tsv"
    params:
        hmms_and_names = config["hmms_and_names"],
        match_directory = output_dir + "/get_matches_viral",
        output_dir = output_dir + "/resolved_matches",
    log:
        logs_dir + "/resolve_conflicts/resolve_conflicts.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/resolve_conflicts.py"