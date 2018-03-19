#!/usr/bin/env python2.7
import argparse
import os
import json
import itertools

import dendropy
from graftm.graftm_package import GraftMPackage

parser = argparse.ArgumentParser()
parser.add_argument('--graftm_package', help='package to look at', required=True)
args = parser.parse_args()

gpkg = GraftMPackage.acquire(args.graftm_package)

taxonomy_hash = gpkg.taxonomy_hash()

taxonomy_to_leaves = {}
for name, taxonomy in taxonomy_hash.items():
    for i in range(len(taxonomy)):
        tax = '; '.join(taxonomy[:(i+1)])
        if tax not in taxonomy_to_leaves:
            taxonomy_to_leaves[tax] = []
        taxonomy_to_leaves[tax].append(name)

refpkg_contents = os.path.join(gpkg.reference_package_path(),'CONTENTS.json')
refpkg = json.loads(open(refpkg_contents).read())
tree_file = os.path.join(gpkg.reference_package_path(),refpkg['files']['tree'])
tree = dendropy.Tree.get(path=tree_file, schema='newick')

print "\t".join([
    'taxonomy',
    'num_lineages',
    'num_lineages_in_mrca',
    'num_disagree_no_taxonomy',
    'euk_count',
    'domain',
    'phylum',
    'class_name',
    'order_name',
    'family',
    'genus',
    'species',
    'eg_bad',
])

label_to_node = {}
for node in tree.leaf_node_iter():
    label_to_node[node.taxon.label] = node

num_printed = 0

for taxonomy, leaves in taxonomy_to_leaves.items():
    if len(leaves) == 1: continue #single member taxonomies are meaningless.
    num_printed += 1
    #if num_printed > 10: break

    leaf_nodes_with_taxonomy = []
    for l in leaves:
        tree_label = l.replace('_',' ')
        if tree_label in label_to_node:
            leaf_nodes_with_taxonomy.append(l)

    if len(leaf_nodes_with_taxonomy) > 0:
        level_to_num_wrong = [0]*8
        mrca = tree.mrca(taxon_labels=list([l.replace('_',' ') for l in leaf_nodes_with_taxonomy]))

        # Calculate how wrong the wrong ones are
        non_belonging_lineages = [l for l in mrca.leaf_iter() if l.taxon.label.replace(' ','_') not in leaf_nodes_with_taxonomy]
        tax_split = taxonomy.split('; ')
        eg_bad = ''
        euks = 0
        for non in non_belonging_lineages:
            num_agreeing_levels = -2
            non_split = taxonomy_hash[non.taxon.label.replace(' ','_')]
            eg_bad = '; '.join(non_split)
            if non_split[1] == 'd__Euryarchaeota':
                euks += 1
            # if eg_bad == 'r__Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Acetobacterales; f__Acetobacteraceae; g__70-18':
            #     import IPython; IPython.embed()
            for i in range(min([len(tax_split),len(non_split)])):
                if tax_split[i] != non_split[i]:
                    num_agreeing_levels = i-2
                    break
            if num_agreeing_levels == len(tax_split) or \
               num_agreeing_levels == len(non_split):
                imcomplete_taxonomies += 1
            num_agreeing_levels += 2
            level_to_num_wrong[num_agreeing_levels] += 1


        print "\t".join(itertools.chain(
            [
                taxonomy,
                str(len(leaf_nodes_with_taxonomy)),
                str(len(mrca.leaf_nodes())),
                str(euks)
            ],
            [str(n) for n in level_to_num_wrong],
            [eg_bad]
        ))


