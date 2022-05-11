import logging
import tempfile
import csv
import numpy as np
import extern

from queue import Queue
from collections import OrderedDict
from .singlem_package import SingleMPackage
from .metapackage import Metapackage

class Condenser:
    """ Combines otu table output for each marker into a single otu table"""

    def condense(self, **kwargs):
        output_otu_table = kwargs.pop('output_otu_table')
        krona_output_file = kwargs.pop('krona')

        to_yield = self.condense_to_otu_table(**kwargs)
        if krona_output_file is not None:
            logging.info("Writing Krona to {}".format(krona_output_file))
            # Cache profiles if we need them twice, otherwise attempt to be
            # low-memory by streaming
            if output_otu_table:
                to_yield = list(to_yield)
            CondensedCommunityProfileKronaWriter.write_krona(to_yield, krona_output_file)
        
        if output_otu_table is not None:
            logging.info("Writing OTU table to {}".format(output_otu_table))
            with open(output_otu_table,'w') as final_otu:
                final_otu.write("\t".join(["sample", "coverage", "taxonomy"])+"\n")
                for condensed_table in to_yield:
                    condensed_table.write_data_to(final_otu)

        logging.info("Finished condense")

    def condense_to_otu_table(self, **kwargs):
        input_otu_table = kwargs.pop('input_streaming_otu_table')
        singlem_packages = kwargs.pop('singlem_packages')
        metapackage_path = kwargs.pop('metapackage_path')
        trim_percent = kwargs.pop('trim_percent') / 100
        min_taxon_coverage = kwargs.pop('min_taxon_coverage')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if singlem_packages and metapackage_path:
            raise Exception("Cannot specify both singlem packages and a metapackage")
        if metapackage_path:
            mpkg = Metapackage.acquire(metapackage_path)
            singlem_package_objects = mpkg.singlem_packages
        elif singlem_packages:
            singlem_package_objects = []
            for path in singlem_packages:
                spkg = SingleMPackage.acquire(path)
                logging.debug("Loading SingleM package: {}".format(spkg.graftm_package_basename()))
                singlem_package_objects.append(spkg)
            logging.info("Loaded %i SingleM packages." % len(singlem_package_objects))
        else:
            # Neither were specified, so use the default set of packages
            logging.debug("Using default set of SingleM packages.")
            mpkg = Metapackage()
            singlem_package_objects = mpkg.singlem_packages


        markers = {} # set of markers used to the domains they target
        target_domains = {"Archaea": [], "Bacteria": [], "Eukaryota": []}
        
        for spkg in singlem_package_objects:
            # ensure v3 packages
            if not spkg.version in [3,4]:
                raise Exception("Only works with v3 or v4 singlem packages.")
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
                
        for domain in target_domains:
            if target_domains[domain] in [1, 2]:
                raise Exception("Number of markers for all domains must either be >= 3 or equal to 0. Only {} markers for domain '{}' found".format(target_domains[domain], domain))

        for sample, sample_otus in input_otu_table.each_sample_otus():
            logging.debug("Processing sample {} ..".format(sample))
            yield self._condense_a_sample(sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage)

    def _condense_a_sample(self, sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage):
        # Stage 1: Build a tree of the observed OTU abundance that is 
        # sample -> gene -> WordNode root
        marker_to_taxon_counts = {} # {sampleID:{gene:wordtree}}}
        excluded_markers = set()

        for otu in sample_otus:
            gene = otu.marker
            coverage = otu.coverage
            tax_split = otu.taxonomy_array()
            # print("Adding {}/{} at coverage {}".format(tax_split,gene,coverage))

            if gene not in markers:
                if gene not in excluded_markers:
                    logging.warning("Gene: {} not in SingleM packages, excluding hits from this marker...".format(gene))
                    excluded_markers.add(gene)
                continue
            
            # ensure OTU is assigned to the domain level or higher
            if len(tax_split) < 2:  # contains at least Root; d__DOMAIN
                continue
            domain = tax_split[1].strip('d__')

            # ensure the domain of this OTU is targeted by the gene
            if domain not in markers[gene]: 
                continue

            if gene not in marker_to_taxon_counts:
                # create new tree for this marker
                marker_to_taxon_counts[gene] = WordNode(None, "Root")

            marker_to_taxon_counts[gene].add_words(tax_split, coverage)

        # Stage 2: Summarise the abundance across the markers for each lineage
        sample_summary_root_node = WordNode(None, "Root")

        for domain, targetted_genes in target_domains.items(): # for each domain
            total_num_markers = len(targetted_genes)

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
                # import IPython; IPython.embed()
                # print()
                # print(node_list)
                # logging.debug("Found abundance {} for taxon {}".format(abundance, node_list[0].get_taxonomy()))
                # if node_list[0].get_taxonomy() in [
                #     ['Root', 'd__Archaea', 'p__Thermoplasmatota', 'c__Poseidoniia'],
                #     ['Root', 'd__Archaea', 'p__Thermoplasmatota', 'c__Poseidoniia', 'o__Poseidoniales'],
                #     ['Root', 'd__Archaea', 'p__Thermoplasmatota', 'c__Poseidoniia', 'o__MGIII']
                # ]:
                #     # 09/03/2021 12:49:57 PM DEBUG: Found abundance 96.214375 for taxon ['Root', 'd__Archaea', 'p__Thermoplasmatota', 'c__Poseidoniia', 'o__MGIII']
                #     [m.get_full_coverage() for m in node_list]
                #     import IPython; IPython.embed()

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

        del marker_to_taxon_counts

        # Stage 3: Correct the coverages by accounting for each node's children
        for node in sample_summary_root_node:
            children_coverage = sum([c.coverage for c in node.children.values()])
            # print("Found cov {} and child coverage {} for {}".format(node.coverage, children_coverage, node.get_taxonomy()))
            if node.word != 'Root':
                node.coverage = node.coverage - children_coverage

                # Apply a general cutoff, which is somewhat arbitrary, but
                # reduces noise. This cutoff also removes the very occasional
                # situations that coverages are negative.
                if node.coverage < min_taxon_coverage:
                    node.coverage = 0

        return CondensedCommunityProfile(sample, sample_summary_root_node)

    def calculate_abundance(self, coverages, total_num_measures, proportiontocut):
        if proportiontocut == 0:
            # the below method does not handle this case correctly
            return sum(coverages) / total_num_measures
        else:
            return _tmean(coverages+([0]*(total_num_measures-len(coverages))), proportiontocut)

def _tmean(data, proportiontocut):
    """
    Returns trimmed mean of an array.
    """
    cut = int(np.floor(len(data) * proportiontocut))
    if cut == 0:
        return np.mean(data)
    else:
        a = sorted(data)
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

    def calculate_level(self):
        '''Return the number of ancestors of this node, which corresponds to taxonomic levels'''
        c = 0
        p = self.parent
        while p is not None:
            c += 1
            p = p.parent
        return c

class CondensedCommunityProfile:
    def __init__(self, sample, tree):
        self.sample = sample
        self.tree = tree

    def write_data_to(self, output_file_io):
        '''Write data to file - IO object is neither opened no closed.'''
        for node in self.tree:
            if round(node.coverage,2) != 0:
                output_file_io.write("\t".join([
                    self.sample,
                    str(round(node.coverage,2)),
                    '; '.join(node.get_taxonomy())
                ])+"\n")

    def breadth_first_iter(self):
        '''Yield a WordNode for each node in the tree.'''
        return self.tree.__iter__()

    @staticmethod
    def each_sample_wise(io):
        '''Read a condensed community profile from an IO, yielding a
        CondenseCommunityProfile object once for each sample.'''

        header = io.readline().strip().split("\t")
        if header != ['sample', 'coverage', 'taxonomy']:
            raise Exception("Unexpected format of condensed community profile file. Expected 'sample', 'coverage', 'taxonomy' as headers.")

        reader = csv.reader(io, delimiter="\t")
        current_sample = None
        current_root = WordNode(None, 'Root')
        taxons_to_wordnode = {current_root.word: current_root}

        for row in reader:
            (sample, coverage, taxonomy) = row
            if sample != current_sample:
                if current_sample is not None:
                    yield CondensedCommunityProfile(current_sample, current_root)
                current_sample = sample
                current_root = WordNode(None, 'Root')
                taxons_to_wordnode = {current_root.word: current_root}
            
            taxons_split = taxonomy.split('; ')
            last_taxon = current_root
            wn = None
            for tax in taxons_split:
                if tax not in taxons_to_wordnode:
                    wn = WordNode(last_taxon, tax)
                    # logging.debug("Adding tax %s with prev %s", tax, last_taxon.word)
                    last_taxon.children[tax] = wn
                    taxons_to_wordnode[tax] = wn
                last_taxon = taxons_to_wordnode[tax]
            wn.coverage = float(coverage)

        if current_sample is not None:
            yield CondensedCommunityProfile(current_sample, current_root)

class CondensedCommunityProfileKronaWriter:
    @staticmethod
    def write_krona(condensed_profiles, output_file):
        cmd = 'ktImportText -o %s' % output_file
        sample_tempfiles = []
        for prof in condensed_profiles:
            f = tempfile.NamedTemporaryFile(prefix='singlem_condense_for_krona',mode='w')
            sample_tempfiles.append(f)

            for node in prof.tree:
                if node.coverage > 0:
                    f.write("\t".join([
                        str(node.coverage),
                        '\t'.join(node.get_taxonomy())
                    ])+"\n")
            f.flush()
            cmd += " %s,'%s'" % (f.name, prof.sample)
        extern.run(cmd)
        for f in sample_tempfiles:
            f.close()




