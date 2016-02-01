import tempfile
import extern
import string
from uc_file import UCFile
from otu_table_entry import OtuTableEntry
from itertools import chain

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
            clustering fraction identity to pass to vsearch e.g. 0.967
            
        Returns
        -------
        None.
        
        Yields
        ------
        Iterates over ClusteredOtu instances, one per sample per cluster,
        in arbitrary order.
        '''
        otus = [o for o in otu_table_collection]
        
        def name_to_index(name):
            return int(string.split(name, ';')[0])

        # write out fasta file numbered to corresponding to the OTU info
        with tempfile.NamedTemporaryFile(prefix='singlem_for_cluster') as f:
            for i, u in enumerate(otus):
                f.write(">%i;size=%i\n" % (i, u.count))
                f.write(u.sequence.replace('-','')+"\n")
            f.flush()
            
            with tempfile.NamedTemporaryFile(prefix='singlem_uc') as uc:
                command = "vsearch --sizein --cluster_size %s --uc %s --id %s" % (f.name, uc.name, str(cluster_identity))
                extern.run(command)
                
                cluster_name_to_sample_to_otus = {}
                for unit in UCFile(open(uc.name)):
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
        