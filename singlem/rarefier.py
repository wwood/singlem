import logging
import copy
import random

from otu_table import OtuTable

class Rarefier:
    def rarefy(self, otu_table_collection, num_to_sample, random_generator=random):
        '''Return an OtuTable rarefied so that only num_to_sample sequences
        are present in each sample. Samples not containing sufficient
        sequences are ignored with a warning.

        This is not a true rarefaction technique because sequences not
        chosen in the rarefaction can still influence the output table
        through the LCA or arbitrary choice operation that has been
        carried out on the input table.

        Also, the rarefier operates on counts rather than predicted
        coverage, skeweing the results toward OTUs that lack
        inserts. But not by a lot, presumably.

        otu_table_collection: OtuTableCollection
            OTU tables iterable
        num_to_sample: int
            number of sequences to sample from each
        '''

        sample_to_gene_to_otu = {}
        to_return = OtuTable()
        for otu in otu_table_collection:
            sample_name = otu.sample_name
            gene = otu.marker
            if sample_name not in sample_to_gene_to_otu:
                sample_to_gene_to_otu[sample_name] = {}
            if gene not in sample_to_gene_to_otu[sample_name]:
                sample_to_gene_to_otu[sample_name][gene] = {}
            if otu.sequence in sample_to_gene_to_otu[sample_name][gene]:
                raise Exception("Found duplicate sequence in OTU table in sample %s, gene %s" % sample_name, gene)
            sample_to_gene_to_otu[sample_name][gene][otu.sequence] = otu

        for sample_name in sample_to_gene_to_otu.keys():
            for gene in sample_to_gene_to_otu[sample_name].keys():
                sequences_to_sample = []
                for sequence, otu in sample_to_gene_to_otu[sample_name][gene].items():
                    for _ in range(otu.count):
                        sequences_to_sample.append(sequence)
                if len(sequences_to_sample) < num_to_sample:
                    logging.warn("Sample %s gene %s only contains %i sequences, so cannot be rarefied. Ignoring this sample/gene combination" % (sample_name, gene, len(sequences_to_sample)))
                    continue
                else:
                    sequences_sampled = random_generator.sample(sequences_to_sample, num_to_sample)
                    sequence_counts = {}
                for seq in sequences_sampled:
                    try:
                        sequence_counts[seq] += 1
                    except KeyError:
                        sequence_counts[seq] = 1

                for seq, count in sequence_counts.items():
                    otu = sample_to_gene_to_otu[sample_name][gene][seq]
                    e = copy.copy(otu)
                    e.count = count
                    to_return.add([e])
        return to_return
