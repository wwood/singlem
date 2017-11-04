import logging
import sys
import math
import squarify
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt


class Appraisal:
    appraisal_results = None

    def plot(self, **kwargs):
        '''Plot this appraisal result to give an idea of what OTUs are accounted for,
        etc.

        Parameters
        ----------
        kwargs:
            output_svg: Filename of SVG to save it to.

        Returns
        -------
        None.
        '''
        output_svg = kwargs.pop('output_svg')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        binned_values = [o.count for o in
                         self.appraisal_results[0].binned_otus if
                         o.marker=='4.15.ribosomal_protein_S2_rpsB']
        binned_values.sort(reverse=True)
        assembled_values = [o.count for o in
                            self.appraisal_results[0].assembled_not_binned_otus() if
                            o.marker=='4.15.ribosomal_protein_S2_rpsB']
        assembled_values.sort(reverse=True)
        not_found_values = [o.count for o in
                            self.appraisal_results[0].not_found_otus if
                            o.marker=='4.15.ribosomal_protein_S2_rpsB']
        not_found_values.sort(reverse=True)
        max_count = max([len(binned_values), len(assembled_values), len(not_found_values)])
        colors = [matplotlib.cm.Pastel1(i) for i in range(max_count)]
        max_area = float(max([sum(binned_values),sum(assembled_values),sum(not_found_values)]))
        fig, axes = plt.subplots(figsize=(12, 10), nrows=3, sharey=True, sharex=True)
        width_and_height = math.sqrt(sum(binned_values)/max_area*100)
        offset = (10 - width_and_height)/2
        axes[0].set_ylim(-offset,10-offset)
        axes[0].set_xlim(-offset,10-offset)
        squarify.plot(binned_values, color=colors, norm_x=width_and_height, norm_y=width_and_height, ax=axes[0])
        width_and_height = math.sqrt(sum(assembled_values)/max_area*100)
        offset = (10 - width_and_height)/2
        axes[1].set_ylim(-offset,10-offset)
        axes[1].set_xlim(-offset,10-offset)
        squarify.plot(assembled_values, color=colors, norm_x=width_and_height, norm_y=width_and_height, ax=axes[1])
        width_and_height = math.sqrt(sum(not_found_values)/max_area*100)
        offset = (10 - width_and_height)/2
        axes[2].set_ylim(-offset,10-offset)
        axes[2].set_xlim(-offset,10-offset)
        squarify.plot(not_found_values, color=colors, norm_x=width_and_height, norm_y=width_and_height, ax=axes[2])
        for a in axes:
            a.set_aspect('equal')
            a.set_ylim(0,10)
            a.set_xlim(0,10)
        #fig.show()

        fig.savefig(output_svg, format='svg')
        #import IPython; IPython.embed()


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
