import logging
import sys
import math
import squarify
import matplotlib
matplotlib.use('Agg') # Must be run the first time matplotlib is imported.
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
            doing_assembly: True or False
            doing_binning: True or False
        Returns
        -------
        None.
        '''
        output_svg_base = kwargs.pop('output_svg_base')
        cluster_identity = kwargs.pop('cluster_identity')
        doing_assembly = kwargs.pop('doing_assembly')
        doing_binning = kwargs.pop('doing_binning')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

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
        for gene in genes:
            self._plot_gene(
                "%s%s.svg" % (output_svg_base, gene),
                cluster_identity,
                gene,
                doing_assembly,
                doing_binning)

    def _plot_gene(self, output_svg, cluster_identity, gene, doing_assembly, doing_binning):
        # Cluster OTUs, making a dict of sequence to cluster object
        plot_info = AppraisalPlotInfo(self, cluster_identity, gene)

        num_samples = len(self.appraisal_results)
        if num_samples == 0:
            raise Exception("Cannot plot an appraisal when there are no samples to appraise")

        fig = plt.figure(figsize=(4.0/5*9 + num_samples*1.0/5*9,5))
        fig.suptitle("SingleM appraisal plot")
        num_panels = 1
        if doing_assembly: num_panels += 1
        if doing_binning: num_panels += 1
        gs = gridspec.GridSpec(num_panels, 4+num_samples)

        legend_axis = fig.add_subplot(gs[:-1,num_samples:])
        self._plot_legend(legend_axis, plot_info)

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
                self._plot_otu(assembled_axis, assembled_values, assembled_colours, max_count,
                               'Unbinned' if sample_number==0 else None)
            unassembled_title = 'Reads' if sample_number==0 else None
            if unassembled_title:
                if unassembled_title and doing_assembly:
                    unassembled_title = 'Unassembled'
                elif doing_binning:
                    unassembled_title = 'Unrecovered'
            self._plot_otu(reads_axis, not_found_values, not_found_colours, max_count,
                           unassembled_title)

            title = sample_appraisal.metagenome_sample_name
            if doing_binning:
                binning_axis.set_title(title)
            elif doing_assembly:
                assembled_axis.set_title(title)
            else:
                reads_axis.set_title(title)

        # Plot the area guide part of the legend
        axis = fig.add_subplot(gs[-1,num_samples])
        self._plot_scale(axis, max_total_count)

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

        axis.bar(x, dy, width=dx, bottom=y, color=colors, align='edge', edgecolor='black', linewidth=0.5)

        axis.set_aspect('equal')
        axis.set_ylim(-0.3,10.3) # Add 0.1 so the edges are not truncated.
        axis.set_xlim(-0.3,10.3)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_axis_off()

        if name is not None:
            fp = matplotlib.font_manager.FontProperties(size=9,weight='bold')
            axis.text(-1, 5,
                      name,
                      font_properties=fp,
                      rotation='vertical',
                      horizontalalignment='center',
                      verticalalignment='center')

    def _plot_legend(self, axis, plot_info):
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
        mycolours = matplotlib.cm.Dark2
        last_index=len(mycolours.colors)
        box_height = 0.6
        space=0.2
        for i, sequence in enumerate(plot_info.ordered_sequences()):
            if i >= last_index: break
            bottom = top-next_y_offset-box_height
            axis.bar(0.1, bottom=bottom, height=box_height, width=0.5,
                     color=mycolours(i), align='edge', edgecolor='black', linewidth=0.5)
            if i==last_index-1:
                t = 'Other'
            else:
                t = plot_info.cluster(sequence).representative_otu.taxonomy.replace('Root; ','')
            axis.text(0.75, top-next_y_offset-box_height*0.6, t, font_properties=fp)
            next_y_offset += box_height + space

    def _plot_scale(self, axis, max_total_count):
        axis.set_aspect('equal')
        axis.set_ylim(-0.3,10.3) # Add 0.1 so the edges are not truncated.
        axis.set_xlim(-0.3,10.3)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_axis_off()
        scale_values = [1, 10, 100] #TODO: Pick scale values so that 111 < max_count
        to_normalise = scale_values + [max_total_count-sum(scale_values)]
        normalised = squarify.normalize_sizes(to_normalise, 10, 10)
        offsets = [9.2, 6.8, 0]
        for i, value in enumerate(scale_values):
            side = math.sqrt(float(value)/max_total_count*100)
            axis.bar(5,bottom=offsets[i],height=side,width=side, color='0.7',edgecolor='black',linewidth=0.5)

            # Add text
            axis.text(10, offsets[i]+side/2, value, verticalalignment='center')
        axis.set_title('OTU count')


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
    def __init__(self, appraisal, cluster_identity, marker):
        '''
        appraisal: Appraisal
        cluster_identity: float, as in Clusterer
        marker: str
            the marker being plotted
        '''
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
            colours.append(matplotlib.cm.Dark2(colour_index))
        return colours

    def cluster(self, sequence):
        return self._sequence_to_cluster[sequence]

    def ordered_sequences(self):
        '''Return representative sequences in same order as the colours'''
        for pair in self._sorted_cluster_rep_and_count:
            yield pair[0]
