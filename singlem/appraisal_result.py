import logging
import sys
import math
import squarify
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt

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
            output_svg: Filename of SVG to save it to.
            cluster_identity: Clusters of OTU calculated threshold (float)
        Returns
        -------
        None.
        '''
        output_svg = kwargs.pop('output_svg')
        cluster_identity = kwargs.pop('cluster_identity')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # Collect OTUs for plotting
        binned_otus = [o for o in
                       self.appraisal_results[0].binned_otus if
                       o.marker=='4.15.ribosomal_protein_S2_rpsB']
        assembled_otus = [o for o in
                          self.appraisal_results[0].assembled_not_binned_otus() if
                          o.marker=='4.15.ribosomal_protein_S2_rpsB']
        not_found_otus = [o for o in
                          self.appraisal_results[0].not_found_otus if
                          o.marker=='4.15.ribosomal_protein_S2_rpsB']

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
        fig, axes = plt.subplots(figsize=(12, 10), nrows=3, sharey=True, sharex=True)
        fig.suptitle(self.appraisal_results[0].metagenome_sample_name)

        self._plot_otu(axes[0], binned_values, binned_colors, max_area)
        self._plot_otu(axes[1], assembled_values, assembled_colors, max_area)
        self._plot_otu(axes[2], not_found_values, not_found_colors, max_area)

        for a in axes:
            a.set_aspect('equal')
            a.set_ylim(0,10)
            a.set_xlim(0,10)
            a.set_xticks([])
            a.set_yticks([])
            a.set_axis_off()
        #fig.show()

        fig.savefig(output_svg, format='svg')
        #import IPython; IPython.embed()

    def _plot_otu(self, axis, sizes, colors, max_area):
        width_and_height = math.sqrt(sum(sizes)/max_area*100)
        squarify.plot(sizes, color=colors, norm_x=width_and_height, norm_y=width_and_height, ax=axis)


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
