import logging
import sys
import math
import re
from textwrap import wrap
import squarify
import matplotlib
matplotlib.use('Agg') # Must be run the first time matplotlib is imported.
import matplotlib.cm
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none' # Export text as text, not paths
import matplotlib.gridspec as gridspec
import numpy

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
            doing_assembly: True or False
            doing_binning: True or False
            gene_to_plot: str name of marker to plot as it appears in the OTU table, or None
                to pick a best representative when single_output_svg is given
            single_output_svg: single output SVG path (not compatible with output_svg_base)
        Returns
        -------
        None.
        '''
        output_svg_base = kwargs.pop('output_svg_base', None)
        cluster_identity = kwargs.pop('cluster_identity')
        doing_assembly = kwargs.pop('doing_assembly')
        doing_binning = kwargs.pop('doing_binning')
        target_gene = kwargs.pop('gene_to_plot', None)
        single_output_svg = kwargs.pop('output_svg', None)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if single_output_svg and output_svg_base:
            raise Exception("Programming error")

        genes = set()
        for r in self.appraisal_results:
            if doing_binning:
                for o in r.binned_otus:
                    genes.add(o.marker)
            if doing_assembly:
                for o in r.assembled_not_binned_otus():
                    genes.add(o.marker)
            for o in r.not_found_otus:
                genes.add(o.marker)

        genes_to_plot = []
        if single_output_svg is not None:
            if target_gene is not None:
                if target_gene in genes:
                    genes_to_plot = [target_gene]
                else:
                    raise Exception("No instances of the marker '%s' found "
                                    "in any OTU tables for plotting." % target_gene)
            else:
                genes_to_plot = [self._pick_representative_marker(doing_assembly, doing_binning)]
        else:
            genes_to_plot = genes

        if len(genes_to_plot) == 1:
            logging.info("Generating plot for marker: %s" % genes_to_plot[0])
        else:
            logging.info("Generating plots for markers: %s" % genes_to_plot)
        for gene in genes_to_plot:
            self._plot_gene(
                "%s%s.svg" % (output_svg_base, gene) if single_output_svg is None \
                else single_output_svg,
                cluster_identity,
                gene,
                doing_assembly,
                doing_binning)

    def _pick_representative_marker(self, doing_assembly, doing_binning):
        '''Pick the marker which has least euclidean distance to the mean
        assembled / binned / metagenome fractions across all the samples.

        '''
        # First collect sample names
        markers = set()
        for result in self.appraisal_results:
            if doing_binning:
                for otu in result.binned_otus: markers.add(otu.marker)
            if doing_assembly:
                for otu in result.assembled_not_binned_otus(): markers.add(otu.marker)
            for otu in result.not_found_otus: markers.add(otu.marker)

        # Calculate geometric mean of each marker for each sample
        marker_to_distances = {}
        for result in self.appraisal_results:
            marker_to_binned = {}
            marker_to_assembled = {}
            marker_to_not_found = {}
            if doing_binning:
                for otu in result.binned_otus:
                    try:
                        marker_to_binned[otu.marker] += otu.count
                    except KeyError:
                        marker_to_binned[otu.marker] = otu.count
            if doing_assembly:
                for otu in result.assembled_not_binned_otus():
                    try:
                        marker_to_assembled[otu.marker] += otu.count
                    except KeyError:
                        marker_to_assembled[otu.marker] = otu.count
            for otu in result.not_found_otus:
                try:
                    marker_to_not_found[otu.marker] += otu.count
                except KeyError:
                    marker_to_not_found[otu.marker] = otu.count
            marker_to_fractions = {}
            for marker in markers:
                numbers = []
                if doing_binning:
                    try:
                        numbers.append(marker_to_binned[marker])
                    except KeyError:
                        numbers.append(0)
                if doing_assembly:
                    try:
                        numbers.append(marker_to_assembled[marker])
                    except KeyError:
                        numbers.append(0)
                try:
                    numbers.append(marker_to_not_found[marker])
                except KeyError:
                    numbers.append(0)
                total = float(sum(numbers))
                marker_to_fractions[marker] = list([float(n)/total for n in numbers])
                logging.debug("Got numbers %s for marker %s in sample %s" % (
                    marker_to_fractions[marker], marker, result.metagenome_sample_name))
            ## Rank markers by distance to the trimmed mean of each of the
            ## numbers as a point in some dimensional space where the dimension
            ## depends on doing_binning etc.
            dimension_means = []
            for i in range(len(marker_to_fractions.values()[0])): # for each dimension
                dimension_means.append(
                    self._trimmean([f for f in marker_to_fractions.values()], 10))
            logging.debug("Found trimmed mean of sample %s as %s" % (
                result.metagenome_sample_name, dimension_means))
            dimension_mean_array = numpy.array(dimension_means)
            for marker in markers:
                distance = numpy.linalg.norm(
                    dimension_mean_array-numpy.array(marker_to_fractions[marker]))
                try:
                    marker_to_distances[marker].append(distance)
                except KeyError:
                    marker_to_distances[marker] = [distance]
        # The best marker is the one with the minimum mean (or equivalently
        # sum) of distances.
        logging.debug("Found overall distances %s" % marker_to_distances)
        return min(marker_to_distances,
                   key = lambda marker: sum(marker_to_distances[marker]))

    def _trimmean(self, arr, percent):
        n = len(arr)
        k = int(round(n*(float(percent)/100)/2))
        return numpy.mean(arr[k+1:n-k])


    def _plot_gene(self, output_svg, cluster_identity, gene, doing_assembly, doing_binning):
        # Cluster OTUs, making a dict of sequence to cluster object
        appraisal_colours = AppraisalPlotColourSet(matplotlib.cm.Dark2)
        plot_info = AppraisalPlotInfo(self, cluster_identity, gene, appraisal_colours)

        num_samples = len(self.appraisal_results)
        if num_samples == 0:
            raise Exception("Cannot plot an appraisal when there are no samples to appraise")

        # margins of figure:
        # 4 + num_samples is the number of x panels
        # 2 or 3 is the number of y panels
        # y needs an extra size for title
        panel_side_length = 1.9
        num_y_panels = 1
        if doing_assembly: num_y_panels += 1
        if doing_binning: num_y_panels += 1
        fig = plt.figure(figsize=(
            (4 + num_samples)*panel_side_length,
            num_y_panels*panel_side_length + 0.5))
        fig.suptitle("SingleM appraisal plot (%s)" % gene)
        gs = gridspec.GridSpec(num_y_panels, 4+num_samples)

        legend_axis = fig.add_subplot(gs[-2,-4])
        self._plot_legend(legend_axis, plot_info, appraisal_colours)

        max_count = plot_info.max_count
        max_total_count = 0
        for sample_number in range(num_samples):
            # Create a different number of panels depending on what is being
            # plotted.
            axis_index = 0
            if doing_binning:
                binning_axis = fig.add_subplot(gs[axis_index,sample_number])
                axis_index += 1
            if doing_assembly:
                assembled_axis = fig.add_subplot(gs[axis_index,sample_number])
                axis_index += 1
            reads_axis = fig.add_subplot(gs[axis_index,sample_number])

            sample_appraisal = self.appraisal_results[sample_number]
            if doing_binning:
                binned_otus = plot_info.sort([o for o in sample_appraisal.binned_otus if o.marker == gene])
            if doing_assembly:
                assembled_otus = plot_info.sort([o for o in sample_appraisal.assembled_not_binned_otus() if o.marker == gene])
            not_found_otus = plot_info.sort([o for o in sample_appraisal.not_found_otus if o.marker == gene])

            if doing_binning:
                binned_values = plot_info.values(binned_otus)
                total = sum(binned_values)
                if total > max_total_count: max_total_count = total
            if doing_assembly:
                assembled_values = plot_info.values(assembled_otus)
                total = sum(assembled_values)
                if total > max_total_count: max_total_count = total
            not_found_values = plot_info.values(not_found_otus)
            total = sum(not_found_values)
            if total > max_total_count: max_total_count = total

            if doing_binning:
                binned_colours = plot_info.colours(binned_otus)
            if doing_assembly:
                assembled_colours = plot_info.colours(assembled_otus)
            not_found_colours = plot_info.colours(not_found_otus)

            if doing_binning:
                self._plot_otu(binning_axis, binned_values, binned_colours, max_count,
                               'Binned' if sample_number==0 else None)
            if doing_assembly:
                subtitle = None
                if sample_number == 0:
                    subtitle = 'Unbinned' if doing_binning else 'Assembled'
                self._plot_otu(assembled_axis, assembled_values,
                               assembled_colours, max_count, subtitle)
            unassembled_title = 'Reads' if sample_number==0 else None
            if unassembled_title:
                if unassembled_title and doing_assembly:
                    unassembled_title = 'Unassembled'
                elif doing_binning:
                    unassembled_title = 'Unrecovered'
            self._plot_otu(reads_axis, not_found_values, not_found_colours, max_count,
                           unassembled_title)

            # Set title
            title = sample_appraisal.metagenome_sample_name
            axis_for_title = None
            if doing_binning: axis_for_title = binning_axis
            elif doing_assembly: axis_for_title = assembled_axis
            else: axis_for_title = reads_axis
            # Wrap so that long titles do not overlap so much
            axis_for_title.set_title("\n".join(wrap(title, 15)))

        # Plot the area guide part of the legend
        axis = fig.add_subplot(gs[-1,num_samples])
        self._plot_scale(axis, max_total_count)

        logging.info("Writing output plot to %s" % output_svg)
        fig.savefig(output_svg, format='svg')



    def _plot_otu(self, axis, sizes, colors, max_count, name):
        width_and_height = math.sqrt(float(sum(sizes))/max_count*100)
        offset = (10 - width_and_height)/2

        normed = squarify.normalize_sizes(sizes, width_and_height, width_and_height)
        rects = squarify.squarify(normed, offset, offset , width_and_height, width_and_height)

        x = [rect['x'] for rect in rects]
        y = [rect['y'] for rect in rects]
        dx = [rect['dx'] for rect in rects]
        dy = [rect['dy'] for rect in rects]

        axis.bar(x, dy, width=dx, bottom=y, color=colors, align='edge', edgecolor='black', linewidth=0.2)

        axis.set_aspect('equal')
        axis.set_ylim(-0.3,10.3) # Add 0.1 so the edges are not truncated.
        axis.set_xlim(-0.3,10.3)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_axis_off()

        if name is not None:
            fp = matplotlib.font_manager.FontProperties(size=9)
            axis.text(-1, 5,
                      name,
                      font_properties=fp,
                      rotation='vertical',
                      horizontalalignment='center',
                      verticalalignment='center')

    def _plot_legend(self, axis, plot_info, appraise_colours):
        axis.set_aspect('equal')
        # Setup global legend properties
        xlim=[-0.3,10.3]
        ylim=xlim
        axis.set_xlim(xlim)
        axis.set_ylim(ylim)
        axis.set_axis_off()
        fp = matplotlib.font_manager.FontProperties(size=9)
        axis.text(0.2, 10.6, 'Cluster taxonomy', font_properties=fp)
        fp = matplotlib.font_manager.FontProperties(size=6)

        # Plot each legend member
        # Plot the 5 clusters with the highest counts
        next_y_offset = 0
        top=10.1
        last_index=len(appraise_colours)
        box_height = 0.8
        space=(top - next_y_offset-len(appraise_colours)*box_height)/(last_index-1)
        for i, sequence in enumerate(plot_info.ordered_sequences()):
            if i >= last_index: break
            bottom = top-next_y_offset-box_height
            axis.bar(0.1, bottom=bottom, height=box_height, width=box_height,
                     color=appraise_colours.colour(i), align='edge', edgecolor='black', linewidth=0.2)
            if i==last_index-1:
                t = 'Other'
            else:
                t = plot_info.cluster(sequence).representative_otu.taxonomy.replace('Root; ','')
                t = re.sub(r'.__','',t)
            axis.text(box_height+0.1+0.2+0.1, top-next_y_offset-box_height*0.6-0.2, t, font_properties=fp)
            next_y_offset += box_height + space

    def _plot_scale(self, axis, max_total_count):
        axis.set_aspect('equal')
        axis.set_ylim(-0.3,10.3) # Add 0.1 so the edges are not truncated.
        axis.set_xlim(-0.3,10.3)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_axis_off()
        scale_values = [1]
        while sum(scale_values)+scale_values[0]*10 < max_total_count:
            scale_values = [scale_values[0]*10] + scale_values
        to_normalise = scale_values + [max_total_count-sum(scale_values)]
        normalised = squarify.normalize_sizes(to_normalise, 10, 10)
        ylim = [0.,10.]
        sides = list([math.sqrt(normalised[i]) for i, value in enumerate(scale_values)])
        overlap = (ylim[1]-ylim[0]-sum(sides))/(len(sides)-1)
        next_bottom = ylim[0]
        xoffset = ylim[1]-ylim[0]-sides[0]
        fp = matplotlib.font_manager.FontProperties(size=6)
        for i, value in enumerate(scale_values):
            side = sides[i]
            axis.bar(10-side/2-0.1-xoffset,
                     bottom=next_bottom,
                     height=side,width=side,
                     color='0.7',edgecolor='black',linewidth=0.2)

            # Add text
            if value == 1:
                text = "1 read"
            else:
                text = "%i reads" % value
            axis.text(10.2-xoffset, next_bottom+side/2, text,
                      verticalalignment='center', font_properties=fp)
            # Setup for next loop iter
            next_bottom += overlap + side
        fp = matplotlib.font_manager.FontProperties(size=9)
        axis.text(0.2, 10.6, 'OTU count', font_properties=fp)


class AppraisalResult:
    num_binned = 0 # binned
    num_assembled = 0 # assembled and/or binned
    num_not_found = 0 # neither assembled nor binned
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

    def num_assembled_not_binned(self):
        count = 0
        for otu in self.assembled_not_binned_otus():
            count += otu.count
        return count


class AppraisalPlotInfo:
    def __init__(self, appraisal, cluster_identity, marker, appraisal_colours):
        '''
        appraisal: Appraisal
        cluster_identity: float, as in Clusterer
        marker: str
            the marker being plotted
        '''
        self.appraisal_colours = appraisal_colours
        logging.debug("Generating plot info for %s" % marker)
        # Collect all OTUs from all samples so that they can be processed
        # together.
        all_binned_otus = []
        all_assembled_not_binned_otus = []
        all_not_found_otus = []
        max_count = 0
        # yuck. Sloppy scope in Python, but not in lambdas when I need it..
        def add_to_totality(otus, totality, max_count):
            count = 0
            for otu in otus:
                if otu.marker==marker:
                    totality.append(otu)
                    count += otu.count
            if count > max_count:
                return count
            else:
                return max_count
        for sample_appraisal in appraisal.appraisal_results:
            max_count = add_to_totality(sample_appraisal.binned_otus, all_binned_otus, max_count)
            max_count = add_to_totality(sample_appraisal.assembled_not_binned_otus(),
                                        all_assembled_not_binned_otus, max_count)
            max_count = add_to_totality(sample_appraisal.not_found_otus, all_not_found_otus, max_count)
        logging.debug("Found maximal count of seqs as %i" % max_count)

        sequence_to_cluster = {}
        cluster_rep_and_count = []
        collection = OtuTableCollection()
        collection.otu_table_objects = [
            all_not_found_otus, all_assembled_not_binned_otus, all_binned_otus
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

        self._sequence_to_cluster = sequence_to_cluster
        self._sorted_cluster_rep_and_count = sorted_cluster_rep_and_count
        self._cluster_sequence_to_order = cluster_sequence_to_order
        self.max_count = max_count

    def sort(self, otus):
        '''Sort OTUs by decreasing count (squarify requires this), then by cluster.
        '''
        return sorted(
            otus,
            key=lambda x: (
                x.count,
                # Use minus so more abundant clusters are first
                -self._cluster_sequence_to_order[self._sequence_to_cluster[x.sequence].sequence]),
            reverse=True)

    def values(self, otus):
        return [otu.count for otu in otus]

    def colours(self, otus):
        '''Return a list of colours in which correspond to the OTUs in the same order
        as self.values().'''
        colours = []
        for otu in otus:
            cluster = self._sequence_to_cluster[otu.sequence]
            colour_index = self._cluster_sequence_to_order[cluster.sequence]
            colours.append(self.appraisal_colours.colour(colour_index))
        return colours

    def cluster(self, sequence):
        return self._sequence_to_cluster[sequence]

    def ordered_sequences(self):
        '''Return representative sequences in same order as the colours'''
        for pair in self._sorted_cluster_rep_and_count:
            yield pair[0]

class AppraisalPlotColourSet:
    '''like colours in matplotlib.cm but have an extra white colour for 'Other'
    OTUs.'''
    def __init__(self, base_colour_scale):
        self.base_colours = base_colour_scale

    def __len__(self):
        return len(self.base_colours.colors)+1

    def colour(self, index):
        if index < len(self.base_colours.colors):
            return self.base_colours(index)
        else:
            return (1.0, 1.0, 1.0, 1.0)
