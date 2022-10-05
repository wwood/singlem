#!/usr/bin/env python3

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2015-2016"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

import argparse
import logging
import sys
import os
import itertools
from graftm.sequence_io import SequenceIO
import queue
import random

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.metapackage import Metapackage

def distance_comparison(singlem_package):
    # Read in taxonomy
    logging.info("Reading taxonomy..")
    gg = singlem_package.graftm_package().taxonomy_hash()
    logging.info("Read in %i taxonomies" % len(gg))
    
    # Read in sequence
    logging.info("Reading sequences..")
    fasta_path = singlem_package.graftm_package().alignment_fasta_path()
    sequences = {seq.name: seq.seq for seq in SequenceIO().read_fasta_file(fasta_path)}
    logging.info("Read in %i sequences" % len(sequences))

    # Remove off-target sequences
    target_domains = singlem_package.target_domains()
    off_target_seq = set()
    for name, taxonomy in gg.items():
        if not taxonomy[0].strip("d__") in target_domains:
            off_target_seq.add(name)
    logging.warn("Found %i off-target sequences, deleting" % len(off_target_seq))
    for name in off_target_seq:
        del gg[name]
    logging.info("After deleting off-target taxonomies now have %i taxes" % len(gg))

    # Create recursive hash of taxonomy ending in sequences
    taxonomic_prefixes = 'd p c o f g s'.split(" ")
    class LineageOrLeaf:
        def __init__(self, parent, name, rank):
            self.parent = parent
            self.children = {} #taxon name to 
            self.name = name
            self.rank = rank
            self.known_to_have_sequence = False
            self.is_leaf = False

        def example_sequence(self):
            # pick an example sequence
            current = self
            while not current.is_leaf:
                current = random.choice([v for v in current.children.values()])
            return random.choice([v for v in current.children.values()])

        def __repr__(self):
            rep = "LineageOrLeage %s with %i children" % (self.name, len(self.children))
            if len(self.children) > 0:
                rep += " e.g. %s" % ([k for k in self.children.keys()][0])
            if self.parent is not None:
                rep += " parent %s" % self.parent.name
            return rep
                                                                  
            
    root = LineageOrLeaf(None, 'root', -1)
    for name, splits in gg.items():
        last = root
        count = 0
        for i, s in enumerate(splits):
            if s[0] != taxonomic_prefixes[i]:
                raise Exception("Didn't expect taxon name %s to not start with '%s'" %(s, taxonomic_prefixes[i]))
            if len(s) < 3: raise Exception("Unexpected taxon %s" % s)
            elif len(s) == 3: break # don't get confused with re-entrant taxonomies
            elif len(s) > 3:
                if s in last.children:
                    last = last.children[s]
                else:
                    current = LineageOrLeaf(last, s, i)
                    last.children[s] = current
                    last = current
            count += 1
        if len(splits) == count and count == len(taxonomic_prefixes):
            last.children[s] = sequences[name]
            last.is_leaf = True
            
    # Prune the LineageOrLeaf tree so that all lineages have at least one sequence
    stack = [root]
    last_rank = len(taxonomic_prefixes)-1
    while len(stack) > 0:
        current = stack.pop()
        if current.known_to_have_sequence:
            continue #children have already been pushed to stack
        if current.rank == last_rank:
            # must have children I think
            if len(current.children) == 0: raise Exception("Leaf nodes should have children")
            else:
                current.known_to_have_sequence = True
                next_parent = current.parent
                while next_parent is not None:
                    next_parent.known_to_have_sequence = True
                    next_parent = next_parent.parent
        elif len(current.children) > 0:
            # push children to stack
            for child in current.children.values():
                stack.append(child)
        else:
            logging.debug("Found orphan taxonomy ending in %s" % current.name)
    # Actually do the pruning. Probably can do this at the same time as above but eh
    stack = [root]
    while len(stack) > 0:
        current = stack.pop()
        if current.is_leaf: continue
        if current.known_to_have_sequence == False:
            logging.debug("Deleting %s" % current.name)
            del current.parent.children[current.name]
        else:
            for child in current.children.values():
                stack.append(child)
    
    
    # Iterate over the taxonomy, calculating distances, calculating a distance between a random choice from each pair of lineages, reporting as we go.
    def compare(seq0, seq1):
        matches = 0
        mismatches = 0
        for i, char in enumerate(seq0):
            char2 = seq1[i]
            if char == '-' and char2 == '':
                mismatches+=1
            elif char == char2:
                matches+=1
            else:
                mismatches+=1
        return float(matches)/(matches+mismatches)
                
    q = queue.Queue()
    q.put(root)
    while not q.empty():
        current = q.get()
        logging.debug("Finding examples from %s" % current)
        
        # choose an example from each pair of lineages
        for pair in itertools.combinations(current.children, 2):
            first = current.children[pair[0]]
            second = current.children[pair[1]]
            example0 = first.example_sequence()
            example1 = second.example_sequence()
            dist = compare(example0, example1)
            if current.rank == -1: break

            print("\t".join([
                singlem_package.graftm_package_basename(),
                taxonomic_prefixes[current.rank],
                current.name,
                str(dist),
                ]))

        if current.rank < last_rank:
            for child in current.children.values():
                q.put(child)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', action="store_true")
    parser.add_argument('--metapackage', help='SingleM metapackage to use', required=True)

    args = parser.parse_args()
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    metapackage = Metapackage.acquire(args.metapackage)

    print("\t".join(["singlem_package", "rank", "taxon", "distance"]))
    for singlem_package in metapackage.singlem_packages:
        distance_comparison(singlem_package)
