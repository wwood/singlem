import logging
import tempfile
import csv
import numpy as np
import extern
import sys
from queue import Queue

from .archive_otu_table import ArchiveOtuTable, ArchiveOtuTableEntry
from .metapackage import Metapackage
from .taxonomy import *

DEFAULT_TRIM_PERCENT = 10
DEFAULT_MIN_TAXON_COVERAGE = 0.35
DEFAULT_GENOME_MIN_TAXON_COVERAGE = 0.1

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
        metapackage_path = kwargs.pop('metapackage_path', None)
        metapackage = kwargs.pop('metapackage', None)

        if metapackage_path:
            logging.info("Using the metapackage at {}".format(metapackage_path))
            metapackage = Metapackage.acquire(metapackage_path)
        elif not metapackage:
            # Neither were specified, so use the default set of packages
            logging.info("Using default SingleM metapackage")
            metapackage = Metapackage.acquire_default()
            
        if metapackage.version < 3:
            raise Exception("Condense function now only works with version 3+ metapackages.")

        # Yield once per sample
        to_yield = self._condense_to_otu_table(metapackage, **kwargs)
        if krona_output_file is not None:
            logging.info("Writing Krona to {}".format(krona_output_file))
            # Cache profiles if we need them twice, otherwise attempt to be
            # low-memory by streaming
            if output_otu_table:
                to_yield = list(to_yield)
            CondensedCommunityProfileKronaWriter.write_krona(to_yield, krona_output_file)
        
        if output_otu_table is not None:
            logging.info("Writing taxonomic profile to {}".format(output_otu_table))
            with open(output_otu_table,'w') as final_otu:
                final_otu.write("\t".join(["sample", "coverage", "taxonomy"])+"\n")
                for condensed_table in to_yield:
                    condensed_table.write_data_to(final_otu)

        logging.info("Finished condense")

    def _condense_to_otu_table(self, metapackage, **kwargs):
        input_otu_table = kwargs.pop('input_streaming_otu_table')
        viral_mode = kwargs.pop('viral_mode', False)
        trim_percent = kwargs.pop('trim_percent', DEFAULT_TRIM_PERCENT) / 100
        min_taxon_coverage = kwargs.pop('min_taxon_coverage', DEFAULT_MIN_TAXON_COVERAGE)
        # apply_expectation_maximisation = kwargs.pop('apply_expectation_maximisation')
        output_after_em_otu_table = kwargs.pop('output_after_em_otu_table', False)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        logging.info("Using minimum taxon coverage of {}".format(min_taxon_coverage))

        markers = {} # set of markers used to the domains they target
        target_domains = {"Archaea": [], "Bacteria": [], "Eukaryota": [], "Viruses": []}
        
        for spkg in metapackage.singlem_packages:
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
                elif domain == "Viruses":
                    target_domains["Viruses"] += [marker_name]
                else:
                    raise Exception("Domain: {} not supported.".format(domain))
                
        for domain in target_domains:
            if target_domains[domain] in [1, 2]:
                raise Exception("Number of markers for all domains must either be >= 3 or equal to 0. Only {} markers for domain '{}' found".format(target_domains[domain], domain))

        for sample, sample_otus in input_otu_table.each_sample_otus(generate_archive_otu_table=True):

            logging.debug("Processing sample {} ..".format(sample))
            apply_diamond_expectation_maximisation = True
            yield self._condense_a_sample(sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage, 
                True, apply_diamond_expectation_maximisation, metapackage, output_after_em_otu_table, viral_mode)

    def _condense_a_sample(self, sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage, 
            apply_query_expectation_maximisation, apply_diamond_expectation_maximisation, metapackage,
            output_after_em_otu_table, viral_mode):


        # Remove off-target OTUs genes
        logging.debug("Total OTU coverage by query: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total OTU coverage by diamond: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))
        logging.info("Removing off-target OTUs from {}".format(sample))
        sample_otus = self._remove_off_target_otus(sample_otus, markers)
        # logging.debug("Total coverage: {}".format(sum([o.coverage for o in sample_otus])))
        logging.info("Total OTU coverage by query: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.info("Total OTU coverage by diamond: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))
        # logging.info("Total coverage: {}".format(sum([o.coverage for o in sample_otus])))

        avg_num_genes_per_species = None
        taxon_marker_counts = None
        if viral_mode:
            logging.info("Running condense in viral mode...")
            if metapackage.version not in [6]:
                raise Exception("Viral mode only works with version 6 metapackages.")
            avg_num_genes_per_species = round(metapackage.avg_num_genes_per_species())
            if avg_num_genes_per_species is None:
                raise Exception("Metapackage does not contain average number of genes per species")
            # taxon_marker_counts = metapackage.get_taxon_marker_counts([";".join(o.taxonomy.split('; ')[1:]) for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])
            query_best_hits = set()
            for o in sample_otus:
                if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD:
                    for best_hit in o.equal_best_hit_taxonomies():
                        query_best_hits.add(best_hit.replace('; ',';')) # pre-emptively strip to avoid issues with whitespace
            # query_best_hits = [o.equal_best_hit_taxonomies() for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD]
            taxon_marker_counts = metapackage.get_taxon_marker_counts(query_best_hits)

        if apply_query_expectation_maximisation:
            sample_otus = self._apply_species_expectation_maximization(sample_otus, trim_percent, target_domains, taxon_marker_counts)
        # logging.info("Total coverage: {}".format(sum([o.coverage for o in sample_otus])))

        if apply_diamond_expectation_maximisation:
            logging.info("Converting DIAMOND IDs to taxons")
            self._convert_diamond_best_hit_ids_to_taxonomies(metapackage, sample_otus)
            # logging.info("Total coverage: {}".format(sum([o.coverage for o in sample_otus])))
            # import pickle; pickle.dump(sample_otus, open("real_data/sample_otus.pkl", "wb"))
            # import pickle; sample_otus = pickle.load(open('real_data/sample_otus.pkl','rb'))
            sample_otus = self._apply_genus_expectation_maximization(sample_otus, target_domains, avg_num_genes_per_species)
            logging.debug("Total OTU coverage: {}".format(sum([o.coverage for o in sample_otus])))

        if output_after_em_otu_table:
            sample_otus.alignment_hmm_sha256s = 'na'
            sample_otus.singlem_package_sha256s = 'na'
            with open(output_after_em_otu_table, 'w') as f:
                sample_otus.write_to(f)

        # Condense via trimmed mean from domain to species
        condensed_otus = self._condense_domain_to_species(sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage, avg_num_genes_per_species, taxon_marker_counts)
        logging.info("Total profile coverage after condense domain to species: {}".format(sum([o.coverage for o in condensed_otus.breadth_first_iter()])))

        # Attribute genus-level down to species level to account for sequencing error
        self._push_down_genus_to_species(condensed_otus, 0.1)
        logging.info("Total profile coverage after push down: {}".format(sum([o.coverage for o in condensed_otus.breadth_first_iter()])))

        self._report_taxonomic_level_assignment_stats(condensed_otus)

        return condensed_otus

    def _report_taxonomic_level_assignment_stats(self, condensed_otus):
        level_coverage = [0.0]*7
        level_count = [0]*7
        for taxon in condensed_otus.breadth_first_iter():
            level = taxon.calculate_level()-1
            if level < 0: continue # skip Root
            level_coverage[level] += taxon.coverage
            if taxon.coverage > 0:
                level_count[level] += 1
        total_coverage = sum(level_coverage)
        logging.info("Taxonomic level coverage:")
        ranks = ['kingdom','phylum','class','order','family','genus','species']
        for level, rank in enumerate(ranks):
            logging.info("{}:\t{:.2f}%\t{} taxons".format(
                rank,
                level_coverage[level]/total_coverage*100 if total_coverage > 0 else 0,
                level_count[level]))

    def _convert_diamond_best_hit_ids_to_taxonomies(self, metapackage, sample_otus):
        '''Input OTU tables assigned taxonomy with diamond have sequence IDs as
        their equal-best hits. This function converts these to taxon strings
        instead.'''

        num_otus_changed = 0
        sequence_ids = set()
        # Step 1: Gather dictionary of sequence IDs to taxon strings
        for otu in sample_otus:
            if otu.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD:
                for seq_id_list in otu.equal_best_hit_taxonomies():
                    for seq_id in seq_id_list:
                        sequence_ids.add(seq_id)

        # Step 2: Get taxon strings
        sequence_id_to_taxon = metapackage.get_taxonomy_of_reads(sequence_ids)

        # Step 3: Write taxons to OTU table
        for otu in sample_otus:
            if otu.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD:
                # Each sequence in the OTU is assigned a separate set of
                # taxon_ids. Maybe we could do something more smart, but for the
                # moment, just assume they are all equally best hits.
                possible_names = set()
                for seq_id_list in otu.equal_best_hit_taxonomies():
                    for seq_id in seq_id_list:
                        taxon_name = sequence_id_to_taxon[seq_id]
                        if not taxon_name[-2].startswith('g__'):
                            if not taxon_name[0] == 'd__Eukaryota':
                                raise Exception("Expected genus level taxon, but found {}, from ID {}".format(taxon_name, seq_id))
                            else:
                                # This can happen when taxonomy is overall
                                # Archaea so not previously filtered out, but
                                # equal-best to Euk
                                logging.debug("Ignoring equal-best hit Eukaryotic taxon {}".format(taxon_name))
                                continue
                        # Record only to genus level
                        if taxon_name[0] != 'Root':
                            taxon_name = ['Root']+taxon_name
                        possible_names.add(';'.join(taxon_name[:-1]))
                otu.data[ArchiveOtuTable.EQUAL_BEST_HIT_TAXONOMIES_INDEX] = list(possible_names)
                num_otus_changed += 1
        logging.info("Converted {} Diamond-assigned OTU taxon_ids to taxon strings".format(num_otus_changed))

    def _push_down_genus_to_species(self, condensed_otus, max_push_down_fraction):
        '''Pushes coverage down proportionally from genus to species level,
        modifying the condensed tree in place'''

        for wordnode in condensed_otus.breadth_first_iter():
            if wordnode.calculate_level() == 6: # Corresponds to genus level
                if not wordnode.word.startswith('g__'):
                    raise Exception("Expected genus level OTU to start with 'g__', found {}".format(wordnode.word))

                species_node_names = list(wordnode.children.keys()) # Do this once to avoid potential changes in order
                total_species_coverage = sum([wordnode.children[s].coverage for s in species_node_names])
                if total_species_coverage == 0:
                    continue # Do nothing, cannot push down at all
                total_genus_coverage = wordnode.coverage + total_species_coverage

                extra_coverage_to_push_down = total_genus_coverage * max_push_down_fraction
                if extra_coverage_to_push_down > wordnode.coverage: # if all genus-only coverage is taken up
                    extra_coverage_to_push_down = wordnode.coverage

                wordnode.coverage -= extra_coverage_to_push_down
                
                for species_node in wordnode.children.values():
                    species_node.coverage += extra_coverage_to_push_down * species_node.coverage / total_species_coverage

    def _condense_domain_to_species(self, sample, sample_otus, markers, target_domains, trim_percent, min_taxon_coverage, avg_num_genes_per_species=None, taxon_marker_counts=None):
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
                taxonomy = node_list[0].get_taxonomy()

                # Calculate stat for this set of markers
                if avg_num_genes_per_species is not None: # Running in viral mode
                    m_coverages = [m.get_full_coverage() for m in node_list]
                    if taxon_marker_counts is not None and len(taxonomy) == 8: # 8 ranks including Root for species-level taxonomy
                        total_num_markers = taxon_marker_counts[';'.join(node_list[0].get_taxonomy())]
                    else:
                        total_num_markers = max(avg_num_genes_per_species, len(m_coverages))
                abundance = self.calculate_abundance(
                    [m.get_full_coverage() for m in node_list],
                    total_num_markers,
                    trim_percent)

                # If stat > 0, add stat to tree and queue children
                if abundance > 0:
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

        # Stage 3: Correct the coverages by accounting for each node's descendants
        for node in sample_summary_root_node:
            children_coverage = sum([c.coverage for c in node.children.values()])
            # print("Found cov {} and child coverage {} for {}".format(node.coverage, children_coverage, node.get_taxonomy()))
            if node.word != 'Root':
                node.coverage = node.coverage - children_coverage

                # Remove the very occasional situations that coverages are
                # negative.
                if node.coverage < 0:
                    node.coverage = 0
                # Apply a general cutoff, which is somewhat arbitrary, but
                # reduces noise. But don't lose that coverage.
                elif node.get_full_coverage() < min_taxon_coverage:
                    # While this node's coverage is set to zero, its parent
                    # should be increased because we don't want to lose coverage.
                    
                    # Each parent node has already been visited, so no need to
                    # calculate the full coverage. Instead, just add this node's
                    # coverage to the most immediate ancestor with non-zero
                    # coverage.
                    parent = node.parent
                    while parent is not None:  # This should never be None since it runs into Root, but just in case
                        if parent.coverage != 0 or parent.word == 'Root':
                            if parent.word != 'Root':
                                # Never add coverage to the Root node
                                parent.coverage += node.coverage
                            break
                        parent = parent.parent
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
        table = ArchiveOtuTable()
        table.fields = sample_otus.fields
        num_no_assignment_otus = sum([otu.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX] for otu in sample_otus if otu.taxonomy_assignment_method() is None])
        num_assigned_otus = sum([otu.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX] for otu in sample_otus if otu.taxonomy_assignment_method() is not None])
        logging.info("Found {:.2f} assigned and {:.2f} unassigned OTU coverage units".format(num_assigned_otus, num_no_assignment_otus))
        if num_no_assignment_otus > num_assigned_otus*0.05:
            logging.warning("Found an expectedly high number of OTUs that have no taxonomy assigned by query or diamond: {} unassigned OTUs and {} assigned, in sample {}".format(num_no_assignment_otus, num_assigned_otus, list(sample_otus)[0].sample_name))
            if num_no_assignment_otus > num_assigned_otus*0.5:
                raise Exception("Stopping: sample {} had too many unassigned OTUs".format(list(sample_otus)[0].sample_name))
        table.data = [otu.data for otu in sample_otus if \
            self._is_targeted_by_marker(otu, otu.taxonomy_array(), markers) and \
            otu.taxonomy_assignment_method() is not None]
        num_no_assignment_otus = sum([otu.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX] for otu in table if otu.taxonomy_assignment_method() is None])
        num_assigned_otus = sum([otu.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX] for otu in table if otu.taxonomy_assignment_method() is not None])
        logging.info("After removing off-target OTUs, found {:.2f} assigned and {:.2f} unassigned OTU coverage units".format(num_assigned_otus, num_no_assignment_otus))
        return table

    def _is_targeted_by_marker(self, otu, tax_split, markers):
        '''return True if the OTU (i.e. domain of the taxonomy) is targeted by
        the marker, else False'''
        if tax_split is None:
            raise Exception("OTU table must contain taxonomy information for condense mode. Was --no-assign-taxonomy specified? It cannot be if the result is to be passed to condense.")
        if tax_split[0] != 'Root':
            raise Exception("OTU tables to condense must contain 'Root' as the first taxon in the taxonomy")
        # ensure OTU is assigned to the domain level or higher
        if len(tax_split) < 2:  # contains at least Root; d__DOMAIN
            return False
        if not tax_split[1].startswith('d__'):
            raise Exception("Unexpected domain name in OTU table: {}".format(tax_split[1]))
        domain = tax_split[1].strip('d__')

        return domain in markers[otu.marker]

    def _apply_genus_expectation_maximization(self, sample_otus, genes_per_domain, avg_num_genes_per_species=None):
        logging.info("Applying genus-wise expectation maximization algorithm to OTU table")
        
        logging.debug("Total coverage by query: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total coverage by diamond: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))
        core_return = self._apply_genus_expectation_maximization_core(sample_otus, 0, genes_per_domain, avg_num_genes_per_species)
        logging.debug("Total coverage by query: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total coverage by diamond: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))

        if core_return is None:
            return sample_otus

        species_to_coverage, best_hit_taxonomy_sets = core_return

        logging.info("Gathering equivalence classes")
        eq_classes = self._gather_equivalence_classes_from_list_of_taxon_lists(best_hit_taxonomy_sets)
        
        logging.info("Demultiplexing OTU table")
        logging.debug("Total coverage by query: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total coverage by diamond: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))
        demux_otus = self._demultiplex_otus(sample_otus, species_to_coverage, eq_classes, DIAMOND_ASSIGNMENT_METHOD)
        logging.debug("Total coverage by query: {}".format(sum([o.coverage for o in demux_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total coverage by diamond: {}".format(sum([o.coverage for o in demux_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))

        logging.info("Finished genus expectation maximization")
        return demux_otus
    
    def _apply_genus_expectation_maximization_core(self, sample_otus, trim_percent, genes_per_domain, avg_num_genes_per_species=None):
        def best_hit_genera_from_otu(otu):
            # The DIAMOND assignments are already truncated to genus level.
            # But the query-based ones go to species level.
            best_hit_taxonomies = otu.equal_best_hit_taxonomies()
            method = otu.taxonomy_assignment_method()
            if method == DIAMOND_ASSIGNMENT_METHOD:
                best_hit_genera = best_hit_taxonomies
            elif method in (QUERY_BASED_ASSIGNMENT_METHOD, QUERY_BASED_ASSIGNMENT_METHOD+'_abandoned'):
                best_hit_genera = set(
                    [';'.join([s.strip() for s in taxon.split(';')[:-1]]) for taxon in best_hit_taxonomies])
            else:
                raise Exception("Unexpected taxonomy assignment method: {}".format(otu.taxonomy_assignment_method()))
            return best_hit_genera

        # Set up initial conditions. The coverage of each genus is set to 1
        genus_to_coverage = {}
        best_hit_taxonomy_sets = set()
        some_em_to_do = False
        for otu in sample_otus:
            best_hit_genera = best_hit_genera_from_otu(otu)
            if otu.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD:
                some_em_to_do = True
            best_hit_taxonomy_sets.add(self._species_list_to_key(best_hit_genera))
            for genus in best_hit_genera:
                if genus not in genus_to_coverage:
                    genus_to_coverage[genus] = 1
        if some_em_to_do is False:
            return None
        logging.debug(best_hit_taxonomy_sets)
        num_steps = 0

        # The fraction of each undecided OTU is the ratio of that class's
        # coverage (coverage in the current iteration) to the total coverage of
        # all best hits of the undecided OTU
        while True: # while not converged
            next_genus_to_gene_to_coverage = {}
            num_steps += 1
            # logging.debug("Starting iteration with species abundances: {}".format(species_to_coverage))

            # Partition out the undecided coverage according to the current
            # iteration's ratios
            for otu in sample_otus:
                unnormalised_coverages = {}
                best_hit_genera = best_hit_genera_from_otu(otu)

                for best_hit_genus in best_hit_genera:
                    if best_hit_genus in genus_to_coverage: # Can get removed during trimmed mean
                        unnormalised_coverages[best_hit_genus] = genus_to_coverage[best_hit_genus]
                total_coverage = sum(unnormalised_coverages.values())
                # if total_coverage == 0:
                #     # This OTU has been removed by the previous iteration through trimmed mean
                #     continue

                for tax, unnormalised_coverage in unnormalised_coverages.items():
                    # Record the total for each gene so a trimmed mean can be taken afterwards
                    if tax not in next_genus_to_gene_to_coverage:
                        next_genus_to_gene_to_coverage[tax] = {}
                    if otu.marker not in next_genus_to_gene_to_coverage[tax]:
                        next_genus_to_gene_to_coverage[tax][otu.marker] = 0
                    next_genus_to_gene_to_coverage[tax][otu.marker] = next_genus_to_gene_to_coverage[tax][otu.marker] + unnormalised_coverage / total_coverage * otu.coverage
                    
            # Calculate the trimmed mean for each genus
            next_genus_to_coverage = {}
            for tax, gene_to_coverage in next_genus_to_gene_to_coverage.items():
                if avg_num_genes_per_species is not None:
                    num_markers = avg_num_genes_per_species
                else:
                    num_markers = len(genes_per_domain[tax.split(';')[1].strip().replace('d__','')])
                logging.debug("Using {} markers for OTU taxonomy {}, with coverages {}".format(num_markers, tax, gene_to_coverage.values()))
                trimmed_mean = self.calculate_abundance(list(gene_to_coverage.values()), num_markers, trim_percent)
                next_genus_to_coverage[tax] = trimmed_mean

            # Has any species changed in abundance by a large enough amount? If not, we're done
            need_another_iteration = False
            for tax, next_coverage in next_genus_to_coverage.items():
                if abs(next_coverage - genus_to_coverage[tax]) > 0.001:
                    logging.debug("Taxonomy {} changed from {} to {}".format(tax, next_coverage, genus_to_coverage[tax]))
                    need_another_iteration = True
                    break
            genus_to_coverage = next_genus_to_coverage
            if not need_another_iteration:
                break
        
        # Round each genome to 4 decimal places in coverage, removing entries with 0 coverage
        # Use 3 decimals to avoid rounding to 0 when one OTU is split between many species
        rounded_genus_to_coverage = {}
        for tax, coverage in genus_to_coverage.items():
            cov2 = round(coverage, 3)
            if cov2 > 0:
                rounded_genus_to_coverage[tax] = cov2
        
        logging.info("Genus-wise EM converged in {} steps".format(num_steps))

        return rounded_genus_to_coverage, \
            list([self._key_to_species_list(k) for k in best_hit_taxonomy_sets])

    def _apply_species_expectation_maximization(self, sample_otus, trim_percent, genes_per_domain, taxon_marker_counts):
        logging.info("Applying species-wise expectation maximization algorithm to OTU table")
        # core_return = self._apply_expectation_maximization_core(sample_otus, trim_percent, genes_per_domain)

        # Do not use trimmed mean for EM, as it seems to give slightly worse
        # results (not well benchmarked though)
        core_return = self._apply_species_expectation_maximization_core(sample_otus, 0, genes_per_domain, taxon_marker_counts=taxon_marker_counts)

        if core_return is None:
            return sample_otus

        species_to_coverage, best_hit_taxonomy_sets = core_return
        logging.debug("Total species_to_coverage coverage: {}".format(sum([cov for cov in species_to_coverage.values()])))
        
        logging.info("Gathering equivalence classes")
        eq_classes = self._gather_equivalence_classes_from_list_of_taxon_lists(best_hit_taxonomy_sets)
        
        logging.debug("Total coverage by query: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total coverage by diamond: {}".format(sum([o.coverage for o in sample_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))
        logging.info("Demultiplexing OTU table")
        
        demux_otus = self._demultiplex_otus(sample_otus, species_to_coverage, eq_classes, QUERY_BASED_ASSIGNMENT_METHOD)
        logging.debug("Total coverage by query: {}".format(sum([o.coverage for o in demux_otus if o.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD])))
        logging.debug("Total coverage by diamond: {}".format(sum([o.coverage for o in demux_otus if o.taxonomy_assignment_method() == DIAMOND_ASSIGNMENT_METHOD])))

        logging.info("Finished expectation maximization")
        return demux_otus

    def _species_list_to_key(self, species_list):
        return '~'.join(species_list)
    def _key_to_species_list(self, key):
        return key.split('~')

    def _apply_species_expectation_maximization_core(self, sample_otus, trim_percent, genes_per_domain, min_genes_for_whitelist=10, proximity_cutoff=0.1, taxon_marker_counts=None):
        # Set up initial conditions. The coverage of each species is set to 1
        species_to_coverage = {}
        best_hit_taxonomy_sets = set()
        some_em_to_do = False
        species_genes = {}


        for otu in sample_otus:
            best_hit_taxonomies = otu.equal_best_hit_taxonomies()
            if otu.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD and best_hit_taxonomies is not None:
                some_em_to_do = True
                best_hit_taxonomy_sets.add(self._species_list_to_key(best_hit_taxonomies))
                for best_hit_tax in best_hit_taxonomies:
                    if best_hit_tax not in species_to_coverage:
                        species_to_coverage[best_hit_tax] = 1
                if len(best_hit_taxonomies) == 1:
                    sp = best_hit_taxonomies[0]
                    if sp not in species_genes:
                        species_genes[sp] = set()
                    species_genes[sp].add(otu.marker)
                species_genes
        if some_em_to_do is False:
            return None
        logging.debug(best_hit_taxonomy_sets)
        species_whitelist = set([sp for (sp, genes) in species_genes.items() if len(genes) >= min_genes_for_whitelist])
        logging.info("Found {} species uniquely hitting >= {} marker genes".format(len(species_whitelist), min_genes_for_whitelist))
        logging.debug("Species whitelist: {}".format(species_whitelist))

        num_steps = 0
        min_num_steps = 50
        # The fraction of each undecided OTU is the ratio of that class's
        # coverage (coverage in the current iteration) to the total coverage of
        # all best hits of the undecided OTU
        while True: # while not converged
            next_species_to_gene_to_coverage = {}
            num_steps += 1
            # logging.debug("Starting iteration with species abundances: {}".format(species_to_coverage))
            for otu in sample_otus:
                unnormalised_coverages = {}
                best_hit_taxonomies = otu.equal_best_hit_taxonomies()
                if otu.taxonomy_assignment_method() == QUERY_BASED_ASSIGNMENT_METHOD and best_hit_taxonomies is not None:
                    for best_hit_tax in best_hit_taxonomies:
                        if best_hit_tax in species_to_coverage: # Can get removed during trimmed mean
                            unnormalised_coverages[best_hit_tax] = species_to_coverage[best_hit_tax]
                total_coverage = sum(unnormalised_coverages.values())
                if total_coverage == 0:
                    # This OTU has been removed by the previous iteration through trimmed mean
                    continue

                for tax, unnormalised_coverage in unnormalised_coverages.items():
                    # Record the total for each gene so a trimmed mean can be taken afterwards
                    if tax not in next_species_to_gene_to_coverage:
                        next_species_to_gene_to_coverage[tax] = {}
                    if otu.marker not in next_species_to_gene_to_coverage[tax]:
                        next_species_to_gene_to_coverage[tax][otu.marker] = 0
                    next_species_to_gene_to_coverage[tax][otu.marker] = next_species_to_gene_to_coverage[tax][otu.marker] + unnormalised_coverage / total_coverage * otu.coverage
                    
            # Calculate the (possibly trimmed) mean for each species
            next_species_to_coverage = {}
            for tax, gene_to_coverage in next_species_to_gene_to_coverage.items():
                if taxon_marker_counts is not None:
                    num_markers = taxon_marker_counts[tax.replace('; ',';')]
                else:
                    num_markers = len(genes_per_domain[tax.split(';')[1].strip().replace('d__','')])
                # logging.debug("Using {} markers for OTU taxonomy {}, with coverages {}".format(num_markers, tax, gene_to_coverage.values()))
                trimmed_mean = self.calculate_abundance(list(gene_to_coverage.values()), num_markers, trim_percent)
                next_species_to_coverage[tax] = trimmed_mean

            # Remove species that appear to be noise based upon having low
            # coverage and proximity to higher coverage species
            if num_steps >= min_num_steps:
                failed_species = self._find_species_with_low_coverage_and_proximity_to_higher_coverage_species(
                    next_species_to_coverage, species_whitelist, proximity_cutoff)
                for failed_s in failed_species:
                    logging.debug("Removing species {} due to low coverage and proximity to higher coverage species".format(failed_species))
                    del next_species_to_coverage[failed_s]

            # Has any species changed in abundance by a large enough amount? If not, we're done
            if num_steps < min_num_steps or len(failed_species) > 0:
                # Always iterate again if we removed any species, because
                # otherwise their coverage contributions will be lost.
                need_another_iteration = True
            else:
                need_another_iteration = False
                for tax, next_coverage in next_species_to_coverage.items():
                    if abs(next_coverage - species_to_coverage[tax]) > 0.001:
                        need_another_iteration = True
                        break

            species_to_coverage = next_species_to_coverage
            if not need_another_iteration:
                break
        
        # Round each genome to 4 decimal places in coverage, removing entries with 0 coverage
        # Use 3 decimals to avoid rounding to 0 when one OTU is split between many species
        rounded_species_to_coverage = {}
        for tax, coverage in species_to_coverage.items():
            cov2 = round(coverage, 3)
            if cov2 > 0:
                rounded_species_to_coverage[tax] = cov2

        logging.info("Species-wise EM converged in {} steps".format(num_steps))

        return rounded_species_to_coverage, \
            list([self._key_to_species_list(k) for k in best_hit_taxonomy_sets])

    def _find_species_with_low_coverage_and_proximity_to_higher_coverage_species(self,
        species_to_coverage, species_whitelist, proximity_cutoff):
        '''Return a list of species that is not in the species whitelist, and
        is less than 10% of the total relative abundance of the genus
        '''
        genus_to_coverage = {}
        for tax, coverage in species_to_coverage.items():
            genus = tax.split(';')[6].strip()
            if genus not in genus_to_coverage:
                genus_to_coverage[genus] = 0
            genus_to_coverage[genus] += coverage

        failed_species = []
        for tax, coverage in species_to_coverage.items():
            genus = tax.split(';')[6].strip()
            if tax not in species_whitelist and coverage < genus_to_coverage[genus] * proximity_cutoff:
                failed_species.append(tax)
        return failed_species


    def _demultiplex_otus(self, sample_otus, species_to_coverage, eq_classes,
    assignment_method):
        ''' Return a new OTU table where the OTUs have been demultiplexed. This
        table likely contains OTUs which have the same window sequence.

        For species that cannot be differentiated according to the eq class,
        collapse them into an LCA taxonomy.
        
        Only OTUs that have been assigned a taxonomy through the given
        assignment_method are considered, the rest are included in the new OTU
        table unchanged.'''

        # Convert eq_classes into a dict of species to LCA
        species_to_equivalence_class_lca = {}
        for sp, eq_class in eq_classes.items():
            species_to_equivalence_class_lca[sp] = TaxonomyUtils.lca_taxonomy_of_strings(eq_class)
        # logging.debug("Species to LCA: {}".format(species_to_equivalence_class_lca))

        # Generate new OTU table. Has to be an Archive because this method is
        # used twice, once for species EM and once for genus EM.
        new_otu_table = ArchiveOtuTable()
        new_otu_table.fields = sample_otus.fields
        for otu in sample_otus:
            if otu.taxonomy_assignment_method() != assignment_method:
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
                logging.debug("LCA to coverage: {}".format(lca_to_coverage))
                total_coverage = sum(lca_to_coverage.values())

                # If there is no LCA, then the speces were removed by the
                # species-wise EM. We do not want to lose the coverage of these
                # so just set them to a genus level hit and keep the coverage.
                if len(lca_to_coverage) == 0:
                    if len(otu.taxonomy_array()) > 7:
                        otu.data[ArchiveOtuTable.TAXONOMY_FIELD_INDEX] = '; '.join(otu.taxonomy_array()[:7])
                    otu.data[ArchiveOtuTable.TAXONOMY_ASSIGNMENT_METHOD_INDEX] += '_abandoned'
                    new_otu_table.add([otu])
                else:
                    for lca, coverage in lca_to_coverage.items():
                        # This used to be a helpful sanity check, but it can
                        # legitimately happen for diamond-assigned OTUs since the
                        # median taxonomy can be longer than the final (and this is
                        # a strlen check not an array length check atm)
                        # if len(lca) < len(otu.taxonomy):
                        #     logging.error("Somehow EM has made taxonomy less specific than the original: {}".format(otu.taxonomy))
                        new_otu = ArchiveOtuTableEntry()
                        new_otu.data = otu.data.copy()
                        new_otu.data[ArchiveOtuTable.TAXONOMY_FIELD_INDEX] = lca
                        new_otu.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX] = coverage / total_coverage * otu.coverage
                        logging.debug("Adding OTU taxonomy {} with coverage {}".format(lca, new_otu.coverage))
                        new_otu_table.add([new_otu])
        return new_otu_table
        
    def _gather_equivalence_classes_from_list_of_taxon_lists(self, species_sets):
        """ Return a dict of species to set of species that are not disinguished
        by any species set """
        species_to_eq_class = {} # Applies equally to taxons generally
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
    CONDENSED_PROFILE_HEADER = ['sample', 'coverage', 'taxonomy']

    def __init__(self, sample, tree):
        self.sample = sample
        self.tree = tree

    @staticmethod
    def write_header_to(output_file_io):
        '''Write header to file - IO object is neither opened nor closed.'''
        output_file_io.write("\t".join(CondensedCommunityProfile.CONDENSED_PROFILE_HEADER)+"\n")

    def write_data_to(self, output_file_io, num_decimals=2):
        '''Write data to file - IO object is neither opened no closed. If num_decimals is None, then the coverage is written as is, without rounding. Otherwise, the coverage is rounded to the number of decimals specified.'''
        for node in self.tree:
            if (num_decimals is None and node.coverage > 0) or round(node.coverage, num_decimals) != 0:
                output_file_io.write("\t".join([
                    self.sample,
                    str(round(node.coverage, num_decimals)) if num_decimals is not None else str(node.coverage),
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
        if header != CondensedCommunityProfile.CONDENSED_PROFILE_HEADER:
            raise Exception("Unexpected format of condensed community profile file. Expected %s as headers." % CondensedCommunityProfile.CONDENSED_PROFILE_HEADER)

        reader = csv.reader(io, delimiter="\t")
        current_sample = None
        current_root = WordNode(None, 'Root')
        taxons_to_wordnode = {current_root.word: current_root}

        for row in reader:
            (sample, coverage, taxonomy) = row
            if float(coverage) == 0:
                continue
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
            if wn is None:
                if taxons_split[-1] in taxons_to_wordnode:
                    # This happens when the profile has more specific ranks
                    # before less specific.
                    wn = taxons_to_wordnode[taxons_split[-1]]
                else:
                    raise Exception("Unexpected processing of taxon {}".format(taxons_split))
            wn.coverage += float(coverage)

        if current_sample is not None:
            yield CondensedCommunityProfile(current_sample, current_root)

    def taxonomic_level_coverage_table(self, assume_8_levels=False):
        '''Return a pl DataFrame with the coverage and relative abundance of
        each taxonomic level. If there are 7 or 8 levels, then the standard
        [root], domain, phylum, etc. levels are assumed. Returning a polars
        dataframe maybe isn't the most pythonic, and so this might be changed in
        the future. But eh for now.'''
        import polars as pl # Import here to avoid requiring it for sandpiper runs
        name_to_coverage = {}
        for node in self.breadth_first_iter():
            node_level = node.calculate_level()
            if node_level == 0:
                continue
            if node_level not in name_to_coverage:
                name_to_coverage[node_level] = 0.
            name_to_coverage[node_level] += node.coverage
        result = pl.DataFrame({
            'level': list(name_to_coverage.keys()),
            'coverage': list(name_to_coverage.values())
        }).with_columns(pl.lit(self.sample).alias('sample')).with_columns(
            ((pl.col('coverage') / pl.col('coverage').sum()).alias('relative_abundance') * 100).round(2),
        )

        if assume_8_levels or len(result.select(pl.col('level')).group_by('level').count()) in [7, 8]:
            # If there's 7 or 8 (including 0) levels, then assume that this is a regular taxonomy going on.
            levels = ['root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            level_id_to_level_name = {i: levels[i] for i in range(len(levels))}
            result = result.with_columns(
                level=pl.col('level').replace_strict(level_id_to_level_name, return_dtype=pl.Utf8)
            )
        # If we assume 8 levels then ensure that each level is present
        if assume_8_levels:
            for level in ['root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                if level not in result['level']:
                    result = result.vstack(pl.DataFrame({
                        'level': [level],
                        'coverage': [0.],
                        'sample': [self.sample],
                        'relative_abundance': [0.]
                    }))
        return result


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




