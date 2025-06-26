#############
### Setup ###
#############
"""
Steps:
Create and activate base environment (env.yml)
Update config.yaml
Run `snakemake --cores 64 --use-conda --retries 2`
"""

localrules:
    all, acquire_and_concat_hmms, hmmsearch_off_target, concat_fscores, resolve_fscores, roundrobin

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
gtdb_proviruses = config["gtdb_proviruses"]

rule all:
    input:
        output_dir + "/roundrobin.done",
        expand(output_dir + "/regenerate/{spkg}.spkg", spkg=hmms_and_names.index)

rule mfqe_viral:
    input:
        touch = output_dir + "/resolve_conflicts.done"
    output:
        touch = output_dir + "/mfqe_viral.done"
    params:
        protein_filepaths = config["viral_faa_list"],
        match_directory = output_dir + "/resolved_matches",
        output_dir = output_dir + "/mfqe_viral",
        logs_dir = logs_dir + "/mfqe_viral"
    threads: workflow.cores
    resources:
        mem_mb = 16 * 1024,
        runtime = 12 * 60
    log:
        logs_dir + "/mfqe_viral.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/mfqe_all.py"

rule transpose_hmms_viral:
    input:
        touch = output_dir + "/mfqe_viral.done"
    output:
        touch = output_dir + "/transpose_hmms_viral.done",
    params:
        protein_filepaths = config["viral_faa_list"],
        matches_dir = output_dir + "/resolved_matches",
        mfqe_dir = output_dir + "/mfqe_viral",
        output_dir = directory(output_dir + "/hmmseq/viral/"), # output is file for each spkg
        taxfiles = [config["viral_tax"]],
        hmms_and_names = output_dir + "/hmms_and_names_noconflict.tsv",
        logs_dir = logs_dir + "/transpose_hmms_viral"
    threads: workflow.cores
    resources:
        mem_mb = 16 * 1024,
        runtime = 8 * 60
    log:
        logs_dir + "/transpose_hmms_viral.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/transpose_hmms_with_sequences_all.py"

rule concatenate_seqs_and_taxonomies_viral:
    input:
        touch = output_dir + "/transpose_hmms_viral.done"
    output:
        done = output_dir + "/hmmseq_concat/viral/{spkg}.done",
        spkg_seq = output_dir + "/hmmseq_concat/viral/{spkg}.faa",
        spkg_tax = output_dir + "/hmmseq_concat/viral/{spkg}_taxonomy.tsv"
    params:
        hmmseq_dir = output_dir + "/hmmseq/viral/",
        concat_dir = output_dir + "/hmmseq_concat/viral",
    shell:
        "mkdir -p {params.concat_dir} && find {params.hmmseq_dir} |grep -F {wildcards.spkg}.faa |parallel --will-cite -j1 --ungroup cat {{}} > {output.spkg_seq} && find {params.hmmseq_dir} |grep -F {wildcards.spkg}_taxonomy.tsv |parallel --will-cite -j1 --ungroup cat {{}} > {output.spkg_tax}; touch {output.done}"

rule graftm_create:
    input:
        done = output_dir + "/hmmseq_concat/viral/{spkg}.done",
        aa_sequences_file = output_dir + "/hmmseq_concat/viral/{spkg}.faa",
        taxonomy_file = output_dir + "/hmmseq_concat/viral/{spkg}_taxonomy.tsv",
    output:
        gpkg = directory(output_dir + "/gpkgs/{spkg}.gpkg"),
        alignment = output_dir + "/gpkgs/{spkg}.gpkg/{spkg}.gpkg.refpkg/{spkg}_deduplicated_aligned.fasta"
    params:
        hmm_arg = lambda wildcards: hmms_and_names.loc[wildcards.spkg, "filepath"],
    threads: 1
    resources:
        mem_mb = 8 * 1024,
        runtime = 24 * 2 * 60
    log:
        logs_dir + "/graftm/{spkg}.gpkg.log"
    conda:
        "envs/singlem.yml"
    shell:
        "graftM create --min_aligned_percent 0 --force --sequences {input.aa_sequences_file} --hmm {params.hmm_arg} --output {output.gpkg} --threads {threads} --no_tree --taxonomy {input.taxonomy_file} 2> {log}"

##########################
# SingleM window finding #
##########################

rule singlem_seqs:
    input:
        alignment = output_dir + "/gpkgs/{spkg}.gpkg/{spkg}.gpkg.refpkg/{spkg}_deduplicated_aligned.fasta",
    output:
        window_position_file = output_dir + "/singlem_seqs/{spkg}.window_position.txt",
    log:
        output_dir + "/singlem_seqs/{spkg}.singlem_seqs.log"
    resources:
        mem_mb = 1 * 1024,
        runtime = 1 * 60
    conda:
        "envs/singlem.yml"
    shell:
        "singlem seqs --alignment {input.alignment} --alignment-type aa > {output.window_position_file} 2> {log}"

# make_initial_spkgs:
rule singlem_create:
    input:
        gpkg = os.path.join(output_dir, "gpkgs", "{spkg}.gpkg"),
        window_position_file = output_dir + "/singlem_seqs/{spkg}.window_position.txt",
    output:
        spkg = directory(output_dir + "/initial_spkgs/{spkg}.spkg")
    log:
        logs_dir + "/singlem_create/{spkg}.singlem_create.log"
    resources:
        mem_mb = 1 * 1024,
        runtime = 1 * 60
    params: 
        taxonomy_file = output_dir + "/hmmseq_concat/viral/{spkg}_taxonomy.tsv",
        gene_description = lambda wildcards: hmms_and_names.loc[wildcards.spkg, "description"],
        target_domains = 'Viruses',
    conda:
        "envs/singlem.yml"
    shell:
        "singlem create --input-graftm-package {input.gpkg} --output {output.spkg} --target-domains {params.target_domains} --hmm-position `cat {input.window_position_file}` --gene-description \"{params.gene_description}\" --force --input-taxonomy {params.taxonomy_file} 2> {log}"

################################
# Acquire off-target sequences #
################################

rule acquire_and_concat_hmms: 
    # fixes weird singlem regenerate crashing from mfqe sequence count mismatch, acquire new search hmms from gpkgs
    input:
        gpkgs = expand(output_dir + "/gpkgs/{spkg}.gpkg", spkg=hmms_and_names.index),
    output:
        hmm_file = output_dir + "/new_hmms_concat.hmm",
        new_hmms_ids = output_dir + "/new_hmms_ids.tsv",
        done = output_dir + "/acquire_and_concat_hmms.done"
    log:
        logs_dir + "/acquire_and_concat_hmms.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/acquire_and_concat_hmms.py"

rule hmmsearch_off_target:
    input:
        done = output_dir + "/acquire_and_concat_hmms.done",
        hmm = output_dir + "/new_hmms_concat.hmm",
        genome_proteins = config["off_target_faa_list"],
    output:
        script_dir = temp(directory(scripts_dir + "/hmmsearch_off_target")),
        touch = output_dir + "/hmmsearch_off_target.done"
    params:
        output_dir = output_dir + "/hmmsearch_off_target",
    threads: workflow.cores
    log:
        logs_dir + "/hmmsearch_off_target.log"
    conda:
        "envs/hmmsearch.yml"
    script:
        "scripts/hmmsearch.py"

rule get_matches_off_target:
    input:
        touch = output_dir + "/hmmsearch_off_target.done"
    output:
        touch = output_dir + "/get_matches_off_target.done"
    params:
        genomes = config["off_target_faa_list"],
        proviruses = gtdb_proviruses,
        hmms_and_names = output_dir + "/new_hmms_ids.tsv",
        hmmsearch_directory = output_dir + "/hmmsearch_off_target",
        output_dir = output_dir + "/get_matches_off_target",
        logs_dir = logs_dir + "/get_matches_off_target"
    threads: workflow.cores
    resources:
        mem_mb = 16 * 1024,
        runtime = 8 * 60
    log:
        logs_dir + "/get_matches_off_target.log"
    conda:
        "envs/get_matches.yml"
    script:
        "scripts/get_matches_all.py"

rule mfqe_off_target:
    input:
        touch = output_dir + "/get_matches_off_target.done"
    output:
        touch = output_dir + "/mfqe_off_target.done"
    params:
        protein_filepaths = config["off_target_faa_list"],
        match_directory = output_dir + "/get_matches_off_target",
        output_dir = output_dir + "/mfqe_off_target",
        logs_dir = logs_dir + "/mfqe_off_target"
    threads: workflow.cores
    resources:
        mem_mb = 16 * 1024,
        runtime = 8 * 60
    log:
        logs_dir + "/mfqe_off_target.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/mfqe_all.py"

rule transpose_hmms_with_offtarget:
    # make sure to clear output directory if rerunning
    input:
        touch = output_dir + "/mfqe_off_target.done"
    output:
        touch = output_dir + "/transpose_hmms_off_target.done"
    params:
        protein_filepaths = config["off_target_faa_list"],
        matches_dir = output_dir + "/get_matches_off_target",
        mfqe_dir= output_dir + "/mfqe_off_target",
        output_dir = directory(output_dir + "/hmmseq/off_target/"), 
        taxfiles = [config["gtdb_arc_tax"], config["gtdb_bac_tax"]],
        hmms_and_names = output_dir + "/new_hmms_ids.tsv",
        logs_dir = logs_dir + "/transpose_hmms_off_target"
    threads: workflow.cores
    resources:
        mem_mb = 16 * 1024,
        runtime = 24 * 60
    log:
        logs_dir + "/transpose_hmms_off_target.log"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/transpose_hmms_with_sequences_all.py"

rule concatenate_seqs_and_taxonomies_off_target:
    input:
        touch = output_dir + "/transpose_hmms_off_target.done"
    output:
        done = output_dir + "/hmmseq_concat/off_target/{spkg}.done",
        spkg_seq = output_dir + "/hmmseq_concat/off_target/{spkg}.faa",
        spkg_tax = output_dir + "/hmmseq_concat/off_target/{spkg}_taxonomy.tsv"
    params:
        hmmseq_dir = output_dir + "/hmmseq/off_target/",
        concat_dir = output_dir + "/hmmseq_concat/off_target",
    resources:
        mem_mb = 8 * 1024,
        runtime = 4 * 60
    shell:
        "mkdir -p {params.concat_dir} && "
        "find {params.hmmseq_dir} -name {wildcards.spkg}.faa |parallel --will-cite -j1 --ungroup cat {{}} > {output.spkg_seq} && "
        "find {params.hmmseq_dir} -name {wildcards.spkg}_taxonomy.tsv |parallel --will-cite -j1 --ungroup cat {{}} > {output.spkg_tax} && touch {output.done}"

rule off_target_dup_rename:
    input:
        done = output_dir + "/hmmseq_concat/off_target/{spkg}.done",
        seqs = output_dir + "/hmmseq_concat/off_target/{spkg}.faa",
        taxonomy = output_dir + "/hmmseq_concat/off_target/{spkg}_taxonomy.tsv",
    output:
        seqs = output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}.faa",
        taxonomy = output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}_taxonomy.tsv",
        done = output_dir + "/hmmseq_concat/off_target_renamed_dups/{spkg}.done",
    params:
        output_dir = output_dir + "/hmmseq_concat/off_target_renamed_dups",
    conda:
        "envs/singlem.yml"
    resources:
        mem_mb = 8 * 1024,
        runtime = 4 * 60
    log:
        log = output_dir + "/logs/off_target_renamed_dups/{spkg}.log"
    script:
        "scripts/rename_off_target_dups.py"

rule singlem_regenerate:
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


##################
# spkg selection #
##################
def get_deduplicated_aligned_fastas(spkg):
    gpkg = spkg.split('/')[-1].rsplit('.',1)[0]
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
        fscores = expand(output_dir + "/fscore/{spkg}.fscore", spkg=hmms_and_names.index),
        done = expand(output_dir + "/fscore/{spkg}.done", spkg=hmms_and_names.index)
    output:
        fscore_list = output_dir + "/fscore_list.tsv",
        done = output_dir + "/concat_fscores.done"
    conda:
        "envs/singlem.yml"
    script:
        "scripts/concat_fscores.py"
        
rule resolve_fscores:
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
    input:
        spkgs = expand(output_dir + "/regenerate/{spkg}.spkg", spkg=hmms_and_names.index),
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
