import os
import csv
import pathlib 
from Bio import Phylo
import logging

tree_path = snakemake.input.tree
viral_fp = snakemake.input.viral_faa_list
outfile = snakemake.output.fscore
BETA = 1.0

# Load viral sequence IDs
viral_ids = set()
with open(viral_fp) as f:
    for line in f:
        viral_ids.add('.'.join(os.path.basename(line.strip()).split('.')[:-1]))

# F-score function
def fscore(tp, fp, fn, beta):
    if tp == 0:
        return 0
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    if precision + recall == 0:
        return 0
    return (1 + beta**2) * precision * recall / (beta**2 * precision + recall)

# Memory-efficient best F-score function
def get_best_fmeasure(tree, beta):
    total_viral = 0
    total_nonviral = 0
    all_leaves = set()

    for leaf in tree.get_terminals():
        if leaf.name.split('~')[1] in viral_ids:
            total_viral += 1
        else:
            total_nonviral += 1
        all_leaves.add(leaf.name)

    best_fscore = -1
    best_leaves = set()

    for branch in tree.get_nonterminals():
        if branch == tree.root:
            continue

        viral_inner = 0
        nonviral_inner = 0
        leaves = set()

        for leaf in branch.get_terminals():
            if leaf.name.split('~')[1] in viral_ids:
                viral_inner += 1
            else:
                nonviral_inner += 1
            leaves.add(leaf.name)

        viral_outer = total_viral - viral_inner
        nonviral_outer = total_nonviral - nonviral_inner

        fscore_inner = fscore(viral_inner, nonviral_inner, total_viral - viral_inner, beta)
        fscore_outer = fscore(viral_outer, nonviral_outer, total_viral - viral_outer, beta)

        if fscore_inner > best_fscore:
            best_fscore = fscore_inner
            best_leaves = leaves
        if fscore_outer > best_fscore:
            best_fscore = fscore_outer
            best_leaves = all_leaves - leaves
    return best_fscore, [best_leaves]

# Multi-iteration F-score pruning
def get_all_fscores(tree, beta=1.0):
    fmeasures = []
    remaining_leaves = set()
    count_viral = 0
    for leaf in tree.get_terminals():
        if leaf.name.split('~')[1] in viral_ids:
            count_viral += 1
        remaining_leaves.add(leaf.name)
    if count_viral == 0:
        return fmeasures
    if len(remaining_leaves) == count_viral:
        return [1.0, 1.0, 1.0], 1
    iterations = 0
    for i in range(5):
        if len(remaining_leaves) < 4 or count_viral == 0:
            break
        fmeasure, clades = get_best_fmeasure(tree, beta)
        fmeasures.append(fmeasure)
        for leaf in clades[0]:
            tree.prune(leaf)
            remaining_leaves.remove(leaf)
            if leaf.endswith('_viral'):
                count_viral -= 1
        iterations += 1
    return fmeasures, iterations

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

if not pathlib.Path(outfile).parent.exists():
    pathlib.Path(outfile).parent.mkdir(parents=True, exist_ok=True)

if not pathlib.Path(snakemake.log[0]).parent.exists():
    pathlib.Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)

try:
    logging.info("Reading tree: " + tree_path)
    tree = Phylo.read(tree_path, "newick")
except ValueError:
    logging.error("Error reading " + tree_path)
fmeasures, iterations = get_all_fscores(tree, BETA)
gene_name = os.path.basename(tree_path)[:-4]
sum_top3 = sum(sorted(fmeasures, reverse=True)[:3])

with open(outfile, 'w') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["gene", "first_fscore", "sum_best_three", "count_iterations"])
    writer.writerow([gene_name, fmeasures[0], sum_top3, iterations])

with open(snakemake.output.done, 'w') as __:
    pass