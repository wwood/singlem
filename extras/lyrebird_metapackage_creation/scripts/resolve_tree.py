import os
import csv
import argparse
from Bio import Phylo
from io import StringIO
import logging

tree_dir = snakemake.input.trees
viral_fp = snakemake.input.viral_faa_list
fscore_list = snakemake.output.fscore_list
resolved_trees_list = snakemake.output.resolved_trees_list
BETA = 1.0

viral_ids = set()
with open(viral_fp) as f:
    for line in f:
        viral_ids.add('.'.join(os.path.basename(line.strip()).split('.')[:-1]))

def fscore(tp, fp, fn, beta):
    if tp == 0:
        return 0
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    if precision + recall == 0:
        return 0
    return (1 + beta**2) * precision * recall / (beta**2 * precision + recall)

# get best f-measure for a tree without rerooting
def get_best_fmeasure(tree, beta):
    f_measures = {}
    total_viral = 0
    total_nonviral = 0
    all_leaves = set()
    for leaf in tree.get_terminals():
        if leaf.name.split('~')[1] in viral_ids:
            total_viral += 1
        else:
            total_nonviral += 1
        all_leaves.add(leaf.name)
    for branch in tree.get_nonterminals():
        if branch == tree.root:
            continue
        viral_count_inner = 0
        nonviral_count_inner = 0
        leaves = set()
        for leaf in branch.get_terminals():
            if leaf.name.split('~')[1] in viral_ids:
                viral_count_inner += 1
            else:
                nonviral_count_inner += 1
            leaves.add(leaf.name)
        viral_count_outer = total_viral - viral_count_inner
        nonviral_count_outer = total_nonviral - nonviral_count_inner

        fscore_inner = fscore(viral_count_inner, nonviral_count_inner, total_viral - viral_count_inner, beta)
        fscore_outer = fscore(viral_count_outer, nonviral_count_outer, total_viral - viral_count_outer, beta)
        if fscore_inner > fscore_outer:
            if fscore_inner not in f_measures:
                f_measures[fscore_inner] = []
            f_measures[fscore_inner].append(leaves)
        else:
            if fscore_outer not in f_measures:
                f_measures[fscore_outer] = []
            f_measures[fscore_outer].append(all_leaves - leaves)
    return max(f_measures.keys()), f_measures[max(f_measures.keys())]

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
        if len(remaining_leaves) < 4:
            break
        if count_viral == 0:
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

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
logging.info("Processing trees")
with open(resolved_trees_list, "w+") as csvfile1:
    with open(fscore_list, "w+") as csvfile2:
        total_trees = 0
        valid_trees = 0
        writer1 = csv.writer(csvfile1, delimiter='\t')
        writer2 = csv.writer(csvfile2, delimiter='\t')
        writer1.writerow(["gene", "first_fscore", "sum_best_three", "count_iterations"])
        writer2.writerow(["gene", "first_fscore", "sum_best_three", "count_iterations"])
        for filename in tree_dir:
            if filename.endswith(".tre"):
                try:
                    tree = Phylo.read(filename, "newick")
                    total_trees += 1
                except ValueError:
                    logging.error("Error reading " + filename)
                    continue
                fmeasures, iterations = get_all_fscores(tree, BETA)
                if len(fmeasures) == 0:
                    continue
                writer2.writerow([os.path.basename(filename)[:-4], fmeasures[0], sum(sorted(fmeasures, reverse=True)[:3]), iterations])
                if len(fmeasures) == 1:
                    writer1.writerow([os.path.basename(filename)[:-4], fmeasures[0], fmeasures[0] * 3, iterations])
                elif fmeasures[0] >= 0.9:
                    writer1.writerow([os.path.basename(filename)[:-4], fmeasures[0], sum(sorted(fmeasures, reverse=True)[:3]), iterations])
                elif fmeasures[0] >= 0.6 and sum(sorted(fmeasures, reverse=True)[:3]) >= 1.8:
                    writer1.writerow([os.path.basename(filename)[:-4], fmeasures[0], sum(sorted(fmeasures, reverse=True)[:3]), iterations])
                else:
                    continue
                valid_trees += 1
        logging.info(f"Processed {total_trees} trees, of which {valid_trees} passed stringency criteria")

with open(snakemake.output.done, 'w') as _: pass
logging.info("Done")
