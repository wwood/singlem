import tempfile
import extern
import string
import logging
from uc_file import UCFile
from otu_table_entry import OtuTableEntry
from otu_table import OtuTable
from itertools import chain
from StringIO import StringIO

class Clusterer:
    def each_cluster(self, otu_table_collection, cluster_identity):
        '''
        Cluster the OTUs in the table collection by sequence identity, using
        preferring sequences with high abundances over those with low abundance.
        Iterate over ClusteredOtu objects that are the result of this.

        Parameters
        ----------
        otu_table_collection: OtuTableCollection
            OTUs to cluster
        cluster_identity: float
            clustering fraction identity to cluster on e.g. 0.967

        Returns
        -------
        None.

        Yields
        ------
        Iterates over ClusteredOtu instances, one per sample per cluster,
        in arbitrary order.
        '''

        # Sort in descending OTU count order so smafa picks OTUs with large
        # counts as the cluster rep.
        otus = sorted([o for o in otu_table_collection], reverse=True, key=lambda x: x.count)

        def name_to_index(name):
            return int(name)

        divergence = None

        # smafa cluster does not accept data streamed in via /dev/stdin, so we
        # make a temporary file.
        with tempfile.NamedTemporaryFile(prefix='singlem_for_cluster') as f:
            for i, u in enumerate(otus):
                if divergence is None:
                    sequence_length = len(u.sequence)
                    divergence = int((1.0-cluster_identity) * sequence_length)
                elif len(u.sequence) != sequence_length:
                    raise Exception(
                        "Currently, for clustering, OTU tables must only contain OTUs that are all equal in length")
                f.write(">%i\n" % i)
                f.write(u.sequence+"\n")
            f.flush()
            uc_contents = extern.run(
                "smafa cluster -d %i '%s'" % (divergence, f.name))
            logging.debug("Found UC file contents from clustering:\n%s" % uc_contents)

            cluster_name_to_sample_to_otus = {}
            for unit in UCFile(StringIO(uc_contents)):
                if unit.target:
                    centre = unit.target
                else:
                    centre = unit.query
                query_otu = otus[name_to_index(unit.query)]
                sample = query_otu.sample_name
                if centre in cluster_name_to_sample_to_otus:
                    if sample in cluster_name_to_sample_to_otus[centre]:
                        cluster_name_to_sample_to_otus[centre][sample].append(query_otu)
                    else:
                        cluster_name_to_sample_to_otus[centre][sample] = [query_otu]
                else:
                    cluster_name_to_sample_to_otus[centre] = {}
                    cluster_name_to_sample_to_otus[centre][sample] = [query_otu]

            for centre_name, sample_to_otus in cluster_name_to_sample_to_otus.items():
                centre = otus[name_to_index(centre_name)]
                ratio = centre.coverage / centre.count

                for sample, sample_otus in sample_to_otus.items():
                    c2 = SampleWiseClusteredOtu()
                    c2.marker = centre.marker
                    c2.sample_name = sample
                    c2.sequence = centre.sequence
                    c2.count = sum([o.count for o in sample_otus])
                    c2.taxonomy = centre.taxonomy
                    c2.coverage = ratio * c2.count

                    c2.data = [
                        c2.marker,
                        c2.sample_name,
                        c2.sequence,
                        c2.count,
                        c2.coverage,
                        c2.taxonomy
                    ]
                    c2.fields = OtuTable.DEFAULT_OUTPUT_FIELDS

                    c2.otus = sample_otus
                    c2.representative_otu = centre
                    yield c2

    def cluster(self, otu_table_collection, cluster_identity):
        '''As per each_cluster(), except that clusters are returned
        as a list of lists of ClusteredOtu objects, so that each cluster is
        together'''
        rep_sequence_to_otus = {}
        for clustered_otu in self.each_cluster(otu_table_collection, cluster_identity):
            repseq = clustered_otu.representative_otu.sequence
            if repseq not in rep_sequence_to_otus:
                rep_sequence_to_otus[repseq] = []
            rep_sequence_to_otus[repseq].append(clustered_otu)

        clusters = []
        for otu_set in rep_sequence_to_otus.values():
            eg = otu_set[0]
            eg.__class__ = Cluster
            eg.otus = list(chain.from_iterable([sample_wise_cluster.otus for sample_wise_cluster in otu_set]))
            clusters.append(eg)

        return Clusters(clusters)

class Cluster(OtuTableEntry):
    '''Basically the same thing as a SampleWiseClusteredOtu except that the
    otus are from all samples, not just a single sample. A semantic
    difference.'''
    pass

class SampleWiseClusteredOtu(Cluster):
    '''A cluster where all of the OTUs are from a single sample, but the
    representative OTU may not be from that sample.'''
    # all otus here are from the sample
    otus = None

    # may or may not be from the sample that self is from
    representative_otu = None

class Clusters:
    def __init__(self, clusters):
        self.clusters = clusters

    def each_otu(self):
        '''Iterate over all OTUs from each clustered_otus'''
        for clustered_otu in self.clustered_otus:
            for otu in clustered_otu.otus:
                yield otu

    def __iter__(self):
        '''Iterate over clusters'''
        for cluster in self.clusters:
            yield cluster
