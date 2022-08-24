import logging
import tempfile
import csv
import numpy as np
import extern
import sys

from queue import Queue

from singlem.otu_table import OtuTable
from singlem.otu_table_entry import OtuTableEntry
from .singlem_package import SingleMPackage
from .metapackage import Metapackage
from .taxonomy import TaxonomyUtils
from .taxonomy import *

# Set CSV field limit to deal with pipe --output-extras as per
# https://github.com/wwood/singlem/issues/89 following
# https://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072
maxInt = sys.maxsize
while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

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
        apply_expectation_maximisation = kwargs.pop('apply_expectation_maximisation')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if singlem_packages and metapackage_path:
            raise Exception("Cannot specify both singlem packages and a metapackage")
        if metapackage_path:
            mpkg = Metapackage.acquire(metapackage_path)
            singlem_package_objects = mpkg.singlem_packages
        else:
            singlem_package_objects = []
            for path in singlem_packages:
                spkg = SingleMPackage.acquire(path)
                logging.debug("Loading SingleM package: {}".format(spkg.graftm_package_basename()))
                singlem_package_objects.append(spkg)
            logging.info("Loaded %i SingleM packages." % len(singlem_package_objects))

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
            yield self._condense_a_sample(sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage, apply_expectation_maximisation)

    def _condense_a_sample(self, sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage, apply_expectation_maximisation):
        if apply_expectation_maximisation:
            # Remove off-target OTUs genes
            sample_otus = self._remove_off_target_otus(sample_otus, markers)
            sample_otus = self._apply_expectation_maximization(sample_otus, trim_percent, target_domains)

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
            
            if not self._is_targeted_by_marker(otu, tax_split, markers):
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

    def _remove_off_target_otus(self, sample_otus, markers):
        """
        Remove OTUs from the OTU table that are not targeted by that marker.
        """
        return [otu for otu in sample_otus if self._is_targeted_by_marker(otu, otu.taxonomy_array, markers)]

    def _is_targeted_by_marker(self, otu, tax_split, markers):
        '''return True if the OTU (i.e. domain of the taxonomy) is targeted by
        the marker, else False'''
        if tax_split[0] != 'Root':
            raise Exception("OTU tables to condense must contain 'Root' as the first taxon in the taxonomy")
        # ensure OTU is assigned to the domain level or higher
        if len(tax_split) < 2:  # contains at least Root; d__DOMAIN
            return False
        if not tax_split[1].startswith('d__'):
            raise Exception("Unexpected domain name in OTU table: {}".format(tax_split[1]))
        domain = tax_split[1].strip('d__')

        return domain in markers[otu.marker]

    def _apply_expectation_maximization(self, sample_otus, trim_percent, genes_per_domain):
        logging.info("Applying core expectation maximization algorithm to OTU table")
        core_return = self._apply_expectation_maximization_core(sample_otus, trim_percent, genes_per_domain)
        if core_return is None:
            return sample_otus

        species_to_coverage, best_hit_taxonomy_sets = core_return
        
        logging.debug("Gathering equivalence classes")
        eq_classes = self._gather_equivalence_classes_from_list_of_species_lists(best_hit_taxonomy_sets)
        
        logging.debug("Demultiplexing OTU table")
        demux_otus = self._demultiplex_otus(sample_otus, species_to_coverage, eq_classes)

        logging.info("Finished expectation maximization")
        return demux_otus

    def _species_list_to_key(self, species_list):
        return '~'.join(species_list)
    def _key_to_species_list(self, key):
        return key.split('~')

    def _apply_expectation_maximization_core(self, sample_otus, trim_percent, genes_per_domain):
        # Set up initial conditions. The coverage of each species is set to 1
        species_to_coverage = {}
        best_hit_taxonomy_sets = set()
        some_em_to_do = False
        for otu in sample_otus:
            best_hit_taxonomies = otu.equal_best_hit_taxonomies()
            if otu.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD and best_hit_taxonomies is not None:
                some_em_to_do = True
                best_hit_taxonomy_sets.add(self._species_list_to_key(best_hit_taxonomies))
                for best_hit_tax in best_hit_taxonomies:
                    if best_hit_tax not in species_to_coverage:
                        species_to_coverage[best_hit_tax] = 1
        if some_em_to_do is False:
            return None
        logging.debug(best_hit_taxonomy_sets)

        # The fraction of each undecided OTU is the ratio of that class's
        # coverage (coverage in the current iteration) to the total coverage of
        # all best hits of the undecided OTU
        while True: # while not converged
            next_species_to_gene_to_coverage = {}
            # logging.debug("Starting iteration with species abundances: {}".format(species_to_coverage))
            for otu in sample_otus:
                unnormalised_coverages = {}
                best_hit_taxonomies = otu.equal_best_hit_taxonomies()
                if otu.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD and best_hit_taxonomies is not None:
                    for best_hit_tax in best_hit_taxonomies:
                        unnormalised_coverages[best_hit_tax] = species_to_coverage[best_hit_tax]
                total_coverage = sum(unnormalised_coverages.values())

                for tax, unnormalised_coverage in unnormalised_coverages.items():
                    # Record the total for each gene so a trimmed mean can be taken afterwards
                    if tax not in next_species_to_gene_to_coverage:
                        next_species_to_gene_to_coverage[tax] = {}
                    if otu.marker not in next_species_to_gene_to_coverage[tax]:
                        next_species_to_gene_to_coverage[tax][otu.marker] = 0
                    next_species_to_gene_to_coverage[tax][otu.marker] = next_species_to_gene_to_coverage[tax][otu.marker] + unnormalised_coverage / total_coverage * otu.coverage
                    
            # Calculate the trimmed mean for each species
            next_species_to_coverage = {}
            for tax, gene_to_coverage in next_species_to_gene_to_coverage.items():
                num_markers = genes_per_domain[tax.split(';')[0]]
                logging.debug("Using {} markers for OTU taxonomy {}, with coverages {}".format(num_markers, tax, gene_to_coverage.values()))
                trimmed_mean = self.calculate_abundance(list(gene_to_coverage.values()), num_markers, trim_percent)
                next_species_to_coverage[tax] = trimmed_mean

            # Has any species changed in abundance by >0.001? If not, we're done
            need_another_iteration = False
            for tax, next_coverage in next_species_to_coverage.items():
                if abs(next_coverage - species_to_coverage[tax]) > 0.00001:
                    need_another_iteration = True
                    break
            species_to_coverage = next_species_to_coverage
            if not need_another_iteration:
                break
        
        # Round each genome to 4 decimal places in coverage, removing entries with 0 coverage
        # Use 4 decimals to avoid rounding to 0 when one OTU is split between many species
        rounded_species_to_coverage = {}
        for tax, coverage in species_to_coverage.items():
            cov2 = round(coverage, 3)
            if cov2 > 0:
                rounded_species_to_coverage[tax] = cov2

        return rounded_species_to_coverage, \
            list([self._key_to_species_list(k) for k in best_hit_taxonomy_sets])

    def _demultiplex_otus(self, sample_otus, species_to_coverage, eq_classes):
        ''' Return a new OTU table where the OTUs have been demultiplexed. This
        table likely contains OTUs which have the same window sequence.

        For species that cannot be differentiated according to the eq class,
        collapse them into an LCA taxonomy '''

        # Convert eq_classes into a dict of species to LCA
        species_to_equivalence_class_lca = {}
        for sp, eq_class in eq_classes.items():
            species_to_equivalence_class_lca[sp] = TaxonomyUtils.lca_taxonomy_of_strings(eq_class)
        # logging.debug("Species to LCA: {}".format(species_to_equivalence_class_lca))

        # Generate new OTU table
        new_otu_table = OtuTable()
        new_otu_table.fields = sample_otus.fields
        for otu in sample_otus:
            if otu.equal_best_hit_taxonomies() != QUERY_BASED_ASSIGNMENT_METHOD:
                new_otu_table.add([otu])
            else:
                lca_to_coverage = {}
                for tax in otu.equal_best_hit_taxonomies():
                    if tax in species_to_coverage: # If not, then it was removed by EM
                        lca = species_to_equivalence_class_lca[tax]
                        cov = species_to_coverage[tax]
                        if lca not in lca_to_coverage:
                            lca_to_coverage[lca] = cov
                        else:
                            lca_to_coverage[lca] += cov
                # logging.debug("LCA to coverage: {}".format(lca_to_coverage))
                total_coverage = sum(lca_to_coverage.values())
                
                for lca, coverage in lca_to_coverage.items():
                    if len(lca) < len(otu.taxonomy):
                        logging.error("Somehow EM has made taxonomy less specific than the original: {}".format(otu.taxonomy))
                    new_otu = OtuTableEntry()
                    new_otu.marker = otu.marker
                    new_otu.sample_name = otu.sample_name
                    new_otu.sequence = otu.sequence
                    new_otu.count = otu.count
                    new_otu.taxonomy = lca
                    new_otu.coverage = coverage / total_coverage * otu.coverage
                    # logging.debug("Adding OTU taxonomy {} with coverage {}".format(lca, new_otu.coverage))
                    new_otu_table.add([new_otu])
        return new_otu_table
        
    def _gather_equivalence_classes_from_list_of_species_lists(self, species_sets):
        """ Return a dict of species to set of species that are not disinguished
        by any species set """
        species_to_eq_class = {}
        for species_set in species_sets:
            # Find all new species in this set. Together they form a new equivalence class
            new_species_in_set = set()
            old_species_in_set = set()
            for species in species_set:
                if species in species_to_eq_class:
                    old_species_in_set.add(species)
                else:
                    new_species_in_set.add(species)
            if len(new_species_in_set) > 0:
                for sp in new_species_in_set:
                    species_to_eq_class[sp] = new_species_in_set

            if len(old_species_in_set) > 0:
                next_set = old_species_in_set
                more_to_go = True
                while more_to_go:
                    one_prev_eq_class = species_to_eq_class[list(next_set)[0]]
                    # Form at most 3 new equivalence classes, corresponding to
                    # the areas of a 2-circle venn diagram
                    #
                    # We use Python set operations
                    extras = one_prev_eq_class - next_set
                    old_others = next_set - one_prev_eq_class
                    intersection = next_set & one_prev_eq_class
                    # Species in the left set form a new class, if there are any
                    if len(extras) > 0:
                        for sp in extras:
                            species_to_eq_class[sp] = extras
                    # Species in the intersection form a new class, if needed
                    if len(intersection) != len(one_prev_eq_class):
                        for sp in intersection:
                            species_to_eq_class[sp] = intersection
                    # Species in the right set are recorded already in one or more other classes. Iterate these classes
                    if len(old_others) > 0:
                        next_set = old_others
                        more_to_go = True
                    else:
                        more_to_go = False
        return species_to_eq_class

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
            
            taxons_split = list([s.strip() for s in taxonomy.split(';')])
            last_taxon = current_root
            wn = None
            logging.debug("Analysing taxon %s" % taxons_split)
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




