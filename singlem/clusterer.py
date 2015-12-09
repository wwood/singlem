import tempfile
import extern

class Clusterer:
    def cluster(self, otu_table_collection):
        '''
        Cluster the 

        Parameters
        ----------
        '''
        otus = [o for o in otu_table_collection]

        # write out fasta file numbered to corresponding to the OTU info
        with tempfile.NamedTemporaryFile(prefix='singlem_for_cluster') as f:
            for i, u in enumerate(otus):
                f.write(">%i;size=%i\n" % (i, u.count))
                f.write(u.sequence)
            f.flush()
            
            with tempfile.NamedTemporaryFile(prefix='singlem_uc') as uc:
                command = "vsearch --sizein --cluster_size %s --uc %s" % (f.name, uc.name)
                extern.run(command)
                
                import IPython; IPython.embed()
                
                # for each line of the uc file
                # get the OTUs that correspond
                # collapse OTUs from the same sample, yield them back
                # take the taxonomy of the representative
                
                
                