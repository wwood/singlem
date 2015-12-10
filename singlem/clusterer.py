import tempfile
import extern
import string
from uc_file import UCFile
from otu_table_entry import OtuTableEntry

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
                        c2 = ClusteredOtu()
                        c2.marker = centre.marker
                        c2.sample_name = sample
                        c2.sequence = centre.sequence
                        c2.count = sum([o.count for o in sample_otus])
                        c2.taxonomy = centre.taxonomy
                        c2.coverage = ratio * c2.count
                        
                        c2.otus = sample_otus
                        c2.representative_otu = centre
                        yield c2
                
class ClusteredOtu(OtuTableEntry):
    # all otus here are from the sample
    otus = None
    
    # may or may not be from the sample that self is from
    representative_otu = None