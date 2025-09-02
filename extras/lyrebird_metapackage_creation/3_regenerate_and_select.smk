#############
### Setup ###
#############
"""
Steps:
Create and activate base environment (env.yml)
Update config.yaml
Run `snakemake --cores 64 --use-conda --retries 2`
"""
localrules: all, concat_fscores, resolve_fscores, roundrobin, cleanup_folders

import pandas as pd
import os

output_dir = config["output_dir"]
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
logs_dir = output_dir + "/logs"
if not os.path.exists(logs_dir):
    os.mkdir(logs_dir)
scripts_dir = output_dir + "/qsub_scripts"

hmms_and_names = pd.read_csv(output_dir + "/hmms_and_names_noconflict.tsv", sep="\t").set_index("gene", drop=False)


def get_filtered_spkgs(wildcards):
    with open(os.path.join(output_dir, "filtered_spkgs.txt")) as f:
        return [line.strip() for line in f if line.strip()]


rule all:
    input:
        output_dir + "/roundrobin.done",
        expand(output_dir + "/regenerate/{spkg}.spkg", spkg=get_filtered_spkgs),
        output_dir + "/cleanup.done",


rule singlem_regenerate:
    """Add off-target sequences to the initial SingleM packages."""
    input:
        off_target_touch = output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}.done",
        singlem_spkg = output_dir + "/initial_spkgs/{spkg}.spkg",
        seqs = output_dir + "/hmmseq_concat/viral/{spkg}.faa",
        taxonomy = output_dir + "/hmmseq_concat/viral/{spkg}_taxonomy.tsv",
        off_target_seqs = output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}.faa",
        off_target_taxonomy = output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}_taxonomy.tsv",
    output:
        spkg = directory(output_dir + "/regenerate/{spkg}.spkg"),
        done = output_dir + "/regenerate/{spkg}.done",
    params:
        sequence_prefix = "{spkg}~",
    resources:
        mem_mb = lambda wildcards, attempt: 32 * 1024 * (2 ** (attempt - 1)),
        runtime = 24 * 2 * 60
    log:
        log = output_dir + "/logs/regenerate/{spkg}.singlem_regenerate.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/regenerate.py"


rule cleanup_folders:
    """Cleans up temporary folders and files created during the workflow."""
    input:
        off_target_rename = expand(output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}.done", spkg=hmms_and_names.index),
    params:
        mfqe_viral = output_dir + "/mfqe_viral",
        mfqe_off_target = output_dir + "/mfqe_off_target",
        hmmseq_viral = output_dir + "/hmmseq/viral",
        hmmseq_off_target = output_dir + "/hmmseq/off_target",
    output:
        done = touch(output_dir + "/cleanup.done"),
    log:
        log = output_dir + "/logs/cleanup.log"
    shell:
        "rm -rf {params.mfqe_viral} {params.mfqe_off_target} {params.hmmseq_viral} {params.hmmseq_off_target} && "
        "rm -rf {output_dir}/hmmsearch_viral {output_dir}/hmmsearch_off_target && "
        "rm -rf {output_dir}/get_matches_viral {output_dir}/get_matches_off_target && "
        "rm -rf {output_dir}/transpose_hmms_viral {output_dir}/transpose_hmms_off_target && "
        "rm -rf {output_dir}/hmmseq_concat/viral {output_dir}/hmmseq_concat/off_target && "
        "rm -rf {output_dir}/hmmseq_concat/off_target_renamed_dups"


#######################################
# Final selection of SingleM packages #
#######################################


def get_deduplicated_aligned_fastas(spkg):
    gpkg = spkg.split('/')[-1].rsplit('.', 1)[0]
    return os.path.join(spkg, gpkg, f"{gpkg}_final.gpkg.refpkg", f"{gpkg}_final_sequences_deduplicated_aligned.fasta")


# would run chainsaw here but chainsaw removes the deduplicated aligned fasta for some reason
rule run_fasttree_mp:
    input:
        done = output_dir + "/regenerate/{spkg}.done",
    params:
        deduplicated_aligned_fasta = get_deduplicated_aligned_fastas(directory(output_dir + "/regenerate/{spkg}.spkg")),
        outdir = output_dir + "/trees"
    output:
        tree = output_dir + "/trees/{spkg}.tre",
    log:
        log = output_dir + "/logs/trees/{spkg}.log"
    benchmark:
        output_dir + "/benchmarking/trees/{spkg}.tre.benchmark"
    threads: lambda wildcards, attempt: 2 ** (attempt - 1)
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * (2 ** (attempt - 1)),
        runtime = 24 * 2 * 60
    conda:
        "envs/singlem.yml"
    shell:
        "mkdir -p {params.outdir} && "
        "OMP_NUM_THREADS={threads} FastTreeMP < {params.deduplicated_aligned_fasta} > {output.tree} 2> {log}"


rule get_fscore:
    """Calculates the F-score fidelity between viral and off-target sequences."""
    input:
        tree = output_dir + "/trees/{spkg}.tre",
        viral_faa_list = config["viral_faa_list"],
    output:
        fscore = temp(output_dir + "/fscore/{spkg}.fscore"),
        done = temp(output_dir + "/fscore/{spkg}.done")
    log:
        log = output_dir + "/logs/fscore/{spkg}.log"
    resources:
        mem_mb = 1 * 1024,
        runtime = 1 * 60
    conda:
        "envs/singlem.yml"
    script:
        "scripts/get_best_fscore.py"


rule concat_fscores:
    input:
        fscores = lambda wildcards: expand(output_dir + "/fscore/{spkg}.fscore", spkg=get_filtered_spkgs(wildcards)),
        done = lambda wildcards: expand(output_dir + "/fscore/{spkg}.done", spkg=get_filtered_spkgs(wildcards))
    output:
        fscore_list = output_dir + "/fscore_list.tsv",
        done = output_dir + "/concat_fscores.done"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/concat_fscores.py"


rule resolve_fscores:
    """Chooses the best HMMs based on the F-score fidelity."""
    input:
        fscore_list = output_dir + "/fscore_list.tsv",
    output:
        resolved_fscores = output_dir + "/resolved_trees_list.tsv",
        done = output_dir + "/resolved_trees.done"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/resolve_fscores.py"


rule roundrobin:
    """Greedy algorithm to maximize viral coverage."""
    input:
        spkgs = lambda wildcards: expand(output_dir + "/regenerate/{spkg}.spkg", spkg=get_filtered_spkgs(wildcards)),
        match_directory = output_dir + "/resolved_matches",
        resolved_trees_list = output_dir + "/resolved_trees_list.tsv",
        done = output_dir + "/resolved_trees.done"
    params:
        hmms_and_names = output_dir + "/hmms_and_names_noconflict.tsv",
    output:
        hmms_and_names_roundrobin = output_dir + "/hmms_and_names_roundrobin.tsv",
        coverages = output_dir + "/roundrobin_species_coverage.tsv",
        done = output_dir + "/roundrobin.done"
    log:
        log = output_dir + "/logs/chosen_spkgs.log"
    resources:
        mem_mb = 4 * 1024,
        runtime = 1 * 60
    conda:
        "envs/singlem.yml"
    script:
        "scripts/roundrobin_search.py"

