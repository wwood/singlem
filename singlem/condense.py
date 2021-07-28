import logging
import tempfile
import csv
import numpy as np
import extern

from queue import Queue
from collections import OrderedDict
from .singlem_package import SingleMPackage

class Condenser:
    """ Combines otu table output for each marker into a single otu table"""

    def condense(self, **kwargs):
        output_otu_table = kwargs.pop('output_otu_table')
        krona_output_file = kwargs.pop('krona')

        condensed_table = self.condense_to_otu_table(**kwargs)

        if output_otu_table is not None:
            condensed_table.write_to_file(output_otu_table)
        if krona_output_file is not None:
            condensed_table.write_krona(krona_output_file)

        logging.info("Finished condense")

    def condense_to_otu_table(self, **kwargs):
        input_otu_table = kwargs.pop('input_otu_table')
        singlem_packages = kwargs.pop('singlem_packages')
        trim_percent = kwargs.pop('trim_percent') / 100
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        markers = {} # set of markers used to the domains they target
        marker_set = set() # set of markers to check if all markers are included
        sample_to_marker_to_taxon_counts = {} # otu tables may contain more than one sample ~ {sampleID:{gene:wordtree}}}
        taxon_counter = {} # count number of genes for each taxonomy before removal
        queues = {} # PriorityQueues for processing taxonomy levels in root to tip order
        target_domains = {"Archaea": [], "Bacteria": [], "Eukaryota": []}
        
        count = 0
        for path in singlem_packages:
            count += 1
            spkg = SingleMPackage.acquire(path)
            logging.debug("Loading SingleM package: {}".format(spkg.graftm_package_basename()))
            # ensure v3 packages
            if not spkg.version in [3]:
                raise Exception("Only works with v3 packages.")
            marker_name = spkg.graftm_package_basename()
            markers[marker_name] = spkg.target_domains()
            # count number of markers for each domain
            for domain in spkg.target_domains():
                if domain == "Archaea":
                    target_domains["Archaea"] += [marker_name]
                elif domain == "Bacteria":
                    target_domains["Bacteria"] += [marker_name]
                elif domain == "Eukaryota":
                    target_domains["Eukaryota"] += [marker_name]
                else:
                    raise Exception("Domain: {} not supported.".format(domain))
        logging.info("Loaded %i SingleM packages." % count)
                
        for domain in target_domains:
            if target_domains[domain] in [1, 2]:
                raise Exception("Number of markers for all domains must either be >= 3 or equal to 0. Only {} markers for domain '{}' found".format(target_domains[domain], domain))

        # Build a tree of the observed OTU abundance that is 
        # sample -> gene -> WordNode root
        logging.info("Reading input OTU table.")
        with open(input_otu_table) as csvfile:
            reader = csv.DictReader(csvfile, delimiter = "\t")
            for line in reader:
                gene = line["gene"]
                sample = line["sample"]
                taxonomy = line["taxonomy"]
                coverage = float(line["coverage"])

                if gene not in markers:
                    if gene not in excluded_markers:
                        logging.info("Gene: {} not in SingleM packages, excluding hits from this marker...".format(gene))
                        excluded_markers.add(gene)
                    continue

                tax_split = taxonomy.split('; ')
                
                # ensure OTU is assigned to the domain level or higher
                if len(tax_split) < 2:  # contains at least Root; d__DOMAIN
                    continue
                domain = tax_split[1].strip('d__')

                # ensure the domain of this OTU is targeted by the gene
                if domain not in markers[gene]: 
                    continue

                if sample not in sample_to_marker_to_taxon_counts:
                    sample_to_marker_to_taxon_counts[sample] = {}
                if gene not in sample_to_marker_to_taxon_counts[sample]:
                    # create new tree for this marker
                    sample_to_marker_to_taxon_counts[sample][gene] = WordNode(None, "Root")

                sample_to_marker_to_taxon_counts[sample][gene].add_words(tax_split, coverage)

        # Summarise the abundance across the markers for each lineage
        sample_to_summarised_taxon_counts = {}
        for sample, marker_to_taxon_counts in sample_to_marker_to_taxon_counts.items():
            sample_summary_root_node = WordNode(None, "Root")
            sample_to_summarised_taxon_counts[sample] = sample_summary_root_node

            # for each domain
            for domain, targetted_genes in target_domains.items():
                total_num_markers = len(target_domains)

                # Extract an initial set of nodes, the marker-wise trees for
                # this domain
                marker_wise_trees = [
                    node for marker, node in marker_to_taxon_counts.items() \
                        if marker in targetted_genes]
                marker_wise_trees_targetted = []
                for marker_tree in marker_wise_trees:
                    for child_name, child_node in marker_tree.children.items():
                        if child_name == 'd__'+domain:
                            marker_wise_trees_targetted.append(child_node)

                # Skip domains if no markers were found
                if len(marker_wise_trees_targetted) == 0:
                    continue

                # Calculate abundances for each level
                node_list_queue = Queue()
                node_list_queue.put(marker_wise_trees_targetted)
                while not node_list_queue.empty():
                    node_list = node_list_queue.get()

                    # Calculate stat for this set of markers
                    abundance = self.calculate_abundance(
                        [m.get_full_coverage() for m in node_list],
                        total_num_markers,
                        trim_percent)
                    # print()
                    # print(node_list)
                    # print("Found abundance {} for taxon {}".format(abundance, node_list[0].get_taxonomy()))

                    # If stat > 0, add stat to tree and queue children
                    if abundance > 0:
                        taxonomy = node_list[0].get_taxonomy()
                        # print(taxonomy, abundance)
                        sample_summary_root_node.add_words(taxonomy, abundance)

                        child_name_to_node_list = {}
                        for tree in node_list:
                            for child in tree.children.values():
                                taxon = child.word
                                if taxon in child_name_to_node_list:
                                    child_name_to_node_list[taxon].append(child)
                                else:
                                    child_name_to_node_list[taxon] = [child]
                        for child_set in child_name_to_node_list.values():
                            # print("adding {}".format(child_set))
                            node_list_queue.put(child_set)

        # Correct the coverages by accounting for each node's children
        for sample, sample_summary_tree in sample_to_summarised_taxon_counts.items():
            for node in sample_summary_tree:
                children_coverage = sum([c.coverage for c in node.children.values()])
                if node.word != 'Root':
                    node.coverage = node.coverage - children_coverage

        return CondensedCommunityProfile(sample_to_summarised_taxon_counts)

    def calculate_abundance(self, coverages, total_num_measures, proportiontocut):
        return _tmean(coverages+([0]*(total_num_measures-len(coverages))), proportiontocut)

def _tmean(data, proportiontocut):
    """
    Returns trimmed mean of an array. 
    
    Slices at least the first and last values of the array.
    """
    a = sorted(data)
    cut = int(np.floor(len(a) * proportiontocut))
    if cut == 0:
        cut = 1
    return np.mean(a[cut:-cut])

class WordNode:
    def __init__(self, parent, word):
        self.parent = parent # WordNode object
        self.word = word
        self.children = {}
        self.coverage = 0
    
    def add_words(self, word_list, coverage):
        """
        Add coverage to this taxonomy to tree.

        Parameters
        ----------
        word_list : LIST
            list of taxonomy terms to be parsed.
        coverage : FLOAT
            coverage to be inserted at the tail end of word_list.

        Raises
        ------
        Exception
            if the current word in word_list mismatches the current node.
        """
        # first word is current
        if self.word == word_list[0]:
            if len(word_list) > 1:
                if word_list[1] in self.children:
                    self.children[word_list[1]].add_words(word_list[1:], coverage)
                else:
                    # word is not present, add to children
                    self.children[word_list[1]] = WordNode(self, word_list[1])
                    self.children[word_list[1]].add_words(word_list[1:], coverage)
            else:
                # last word
                self.coverage += coverage
        else:
            raise Exception("Word %s does not match current node %s" % (word_list[0]), self.word)
    
    def get_full_coverage(self):
        """
        Get full coverage sum of this node and all its children.
        """
        covs = 0
        covs += self.coverage
        if len(self.children) == 0: # leaf node
            return covs
        for node in self.children.values():
            covs += node.get_full_coverage()
        return covs

    def get_taxonomy(self):
        if self.parent is None:
            return ['Root']
        else:
            return self.parent.get_taxonomy() + [self.word]

    def __iter__(self):
        self.iter_queue = Queue()
        self.iter_queue.put(self)
        return self

    def __next__(self):
        if self.iter_queue.empty():
            raise StopIteration()
        else:
            node = self.iter_queue.get()
            for c in node.children.values():
                self.iter_queue.put(c)
            return node

class CondensedCommunityProfile:
    def __init__(self, sample_to_tree):
        self.sample_to_tree = sample_to_tree

    def write_to_file(self, output_file):
        with open(output_file,'w') as final_otu:
            final_otu.write("\t".join(["sample", "coverage", "taxonomy"])+"\n")
            for sample, sample_summary_tree in self.sample_to_tree.items():
                for node in sample_summary_tree:
                    if node.coverage != 0:
                        final_otu.write("\t".join([
                            sample,
                            str(node.coverage.round(2)),
                            '; '.join(node.get_taxonomy())
                        ])+"\n")

    def write_krona(self, output_file):
        cmd = 'ktImportText -o %s' % output_file
        sample_tempfiles = []
        for sample, tree in self.sample_to_tree.items():
            f = tempfile.NamedTemporaryFile(prefix='singlem_condense_for_krona',mode='w')
            sample_tempfiles.append(f)

            for node in tree:
                if node.coverage > 0:
                    f.write("\t".join([
                        str(node.coverage),
                        '\t'.join(node.get_taxonomy())
                    ])+"\n")
            f.flush()
            cmd += " %s,'%s'" % (f.name, sample)
        extern.run(cmd)
        for f in sample_tempfiles:
            f.close()




