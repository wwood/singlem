import logging
import sys
import math
import squarify
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from clusterer import Clusterer
from otu_table_collection import OtuTableCollection

class Appraisal:
    appraisal_results = None

    def plot(self, **kwargs):
        '''Plot this appraisal result to give an idea of what OTUs are accounted for,
        etc.

        Parameters
        ----------
        kwargs:
            output_svg_base: Basename of SVGs (one SVG per gene is generated)
            cluster_identity: Clusters of OTU calculated threshold (float)
        Returns
        -------
        None.
        '''
        output_svg_base = kwargs.pop('output_svg_base')
        cluster_identity = kwargs.pop('cluster_identity')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        genes = set()
        for r in self.appraisal_results:
            for o in r.binned_otus:
                genes.add(o.marker)
            for o in r.assembled_not_binned_otus():
                genes.add(o.marker)
            for o in r.not_found_otus:
                genes.add(o.marker)
        for gene in genes:
            self._plot_gene(
                "%s%s.svg" % (output_svg_base, gene),
                cluster_identity,
                gene)

    def _plot_gene(self, output_svg, cluster_identity, gene):
        # Collect OTUs for plotting
        binned_otus = [o for o in
                       self.appraisal_results[0].binned_otus if
                       o.marker==gene]
        assembled_otus = [o for o in
                          self.appraisal_results[0].assembled_not_binned_otus() if
                          o.marker==gene]
        not_found_otus = [o for o in
                          self.appraisal_results[0].not_found_otus if
                          o.marker==gene]

        # Cluster OTUs, making a dict of sequence to cluster object
        sequence_to_cluster = {}
        cluster_rep_and_count = []
        collection = OtuTableCollection()
        collection.otu_table_objects = [
            not_found_otus, assembled_otus, binned_otus
        ]
        for cotu in Clusterer().cluster(collection, cluster_identity):
            cluster_rep_and_count.append([cotu.sequence, cotu.count])
            for otu in cotu.otus:
                sequence_to_cluster[otu.sequence] = cotu

        # Sort the OTUs by descending order of counts, so that more abundant
        # OTUs get colour.
        sorted_cluster_rep_and_count = sorted(
            cluster_rep_and_count, key=lambda x: x[1], reverse=True)
        cluster_sequence_to_order = {}
        i = 0
        for pair in sorted_cluster_rep_and_count:
            cluster_sequence_to_order[pair[0]] = i
            i += 1

        # Sort OTUs by decreasing count (squarify requires this), then by
        # cluster.
        def sort_otus(otus):
            return sorted(
                otus,
                key=lambda x: (
                    x.count,
                    # Use minus so more abundant clusters are first
                    -cluster_sequence_to_order[sequence_to_cluster[x.sequence].sequence]),
                reverse=True)
        binned_otus = sort_otus(binned_otus)
        assembled_otus = sort_otus(assembled_otus)
        not_found_otus = sort_otus(not_found_otus)

        # Calculate counts
        binned_values = [o.count for o in binned_otus]
        assembled_values = [o.count for o in assembled_otus]
        not_found_values = [o.count for o in not_found_otus]

        # Calculate colours
        def colours(otus):
            colours = []
            for otu in otus:
                cluster = sequence_to_cluster[otu.sequence]
                colour_index = cluster_sequence_to_order[cluster.sequence]
                colours.append(matplotlib.cm.Pastel1(colour_index))
            return colours
        binned_colors = colours(binned_otus)
        assembled_colors = colours(assembled_otus)
        not_found_colors = colours(not_found_otus)

        max_count = max([len(binned_values), len(assembled_values), len(not_found_values)])
        max_area = float(max([sum(binned_values),sum(assembled_values),sum(not_found_values)]))

        fig = plt.figure(figsize=(9,5))
        fig.suptitle("Appraisal plot for %s" %\
                     self.appraisal_results[0].metagenome_sample_name)
        gs = gridspec.GridSpec(3,5)

        binning_axis = fig.add_subplot(gs[0,0])
        assembled_axis = fig.add_subplot(gs[1,0])
        reads_axis = fig.add_subplot(gs[2,0])
        legend_axis = fig.add_subplot(gs[:,1:])

        self._plot_otu(binning_axis, binned_values, binned_colors, max_area, 'Binned')
        self._plot_otu(assembled_axis, assembled_values, assembled_colors, max_area, 'Unbinned')
        self._plot_otu(reads_axis, not_found_values, not_found_colors, max_area, 'Unassembled')

        self._plot_legend(legend_axis, sorted_cluster_rep_and_count, sequence_to_cluster)

        fig.tight_layout()
        fig.savefig(output_svg, format='svg')


    def _plot_otu(self, axis, sizes, colors, max_area, name):
        width_and_height = math.sqrt(sum(sizes)/max_area*100)
        offset = (10 - width_and_height)/2

        normed = squarify.normalize_sizes(sizes, width_and_height, width_and_height)
        rects = squarify.squarify(normed, offset, offset , width_and_height, width_and_height)

        x = [rect['x'] for rect in rects]
        y = [rect['y'] for rect in rects]
        dx = [rect['dx'] for rect in rects]
        dy = [rect['dy'] for rect in rects]

        axis.bar(x, dy, width=dx, bottom=y, color=colors, align='edge', edgecolor='black')

        axis.set_aspect('equal')
        axis.set_ylim(-0.3,10.3) # Add 0.1 so the edges are not truncated.
        axis.set_xlim(-0.3,10.3)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_axis_off()

        fp = matplotlib.font_manager.FontProperties(size=9,weight='bold')
        axis.text(-1, width_and_height/2,name,font_properties=fp,rotation='vertical',horizontalalignment='center',
                  verticalalignment='center')

    def _plot_legend(self, axis, sorted_cluster_rep_and_count, sequence_to_cluster):
        # Setup global legend properties
        axis.set_xlim([0,10])
        axis.set_ylim([0,14])
        axis.set_axis_off()
        fp = matplotlib.font_manager.FontProperties(size=9,weight='bold')
        axis.text(0.2, 13, 'Taxonomy', font_properties=fp)
        fp = matplotlib.font_manager.FontProperties(size=6)

        # Plot each legend member
        # Plot the 5 clusters with the highest counts
        next_y_offset = 0
        num_to_print = 5
        top=12.6
        last_index=9 # From Pastel1
        box_height = 0.6
        space=0.2
        for i in range(last_index):
            bottom = top-next_y_offset-box_height
            axis.bar(0.1, bottom=bottom, height=box_height, width=0.5,
                     color=matplotlib.cm.Pastel1(i), align='edge', edgecolor='black')
            if i==last_index-1:
                t = 'Other'
            else:
                t = sequence_to_cluster[sorted_cluster_rep_and_count[i][0]].taxonomy.replace('Root; ','')
            axis.text(0.75, top-next_y_offset-box_height*0.6, t, font_properties=fp)
            next_y_offset += box_height + space



class AppraisalResult:
    num_binned = 0
    num_assembled = 0
    num_not_found = 0
    metagenome_sample_name = None

    def __init__(self):
        self.binned_otus = []
        self.assembled_otus = []
        self.not_found_otus = []

    def assembled_not_binned_otus(self):
        '''Iterate (yield) over those OTUs that are assembled but not binned.
        '''
        binned_otu_sequences = set()
        for otu in self.binned_otus:
            binned_otu_sequences.add(otu.sequence)
        for otu in self.assembled_otus:
            if otu.sequence not in binned_otu_sequences:
                yield otu
