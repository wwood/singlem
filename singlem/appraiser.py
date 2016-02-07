import logging
import sys

from clusterer import Clusterer
from otu_table_collection import OtuTableCollection

class Appraiser:
    def appraise(self, **kwargs):
        '''Given a collection of OTU tables derived from samples, and OTU
        table(s) corresponding to a collection of recovered genomes, how
        much of the community has been recovered in those genomes?
        
        Returns
        -------
        An Appraisal object containing appraisals for each metagenome
        '''
        genome_otu_table_collection = kwargs.pop('genome_otu_table_collection')
        metagenome_otu_table_collection = kwargs.pop('metagenome_otu_table_collection')
        cluster_identity = kwargs.pop('cluster_identity', None)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        if cluster_identity is None:
            # read in genome OTU sequences
            genome_otu_sequences = set()
            genome_names = set()
            for otu in genome_otu_table_collection:
                genome_otu_sequences.add(otu.sequence)
                genome_names.add(otu.sample_name)
            logging.info("Read in %i unique sequences from the %i reference genomes" % (len(genome_otu_sequences), len(genome_names)))
            
            # read in metagenome OTU sequences
            sample_name_to_appraisal = {}
            for otu in metagenome_otu_table_collection:
                try:
                    appraisal = sample_name_to_appraisal[otu.sample_name]
                except KeyError:
                    appraisal = AppraisalResult()
                    appraisal.metagenome_sample_name = otu.sample_name
                    sample_name_to_appraisal[otu.sample_name] = appraisal
                    
                count = otu.count
                if otu.sequence in genome_otu_sequences:
                    appraisal.num_found += count
                else:
                    appraisal.num_not_found += count
                    
            app = Appraisal()
            app.appraisal_results = sample_name_to_appraisal.values()
            return app
        
        else:
            # otu collections may be streams, but we need to access them
            # twice, so slurp.
            metagenome_otu_table = OtuTableCollection()
            metagenome_otu_table.otu_table_objects.append(list(metagenome_otu_table_collection))
            logging.info("Read in %i OTU sequences from metagenomes" % len(list(metagenome_otu_table)))
            genomes_otu_table = OtuTableCollection()
            genomes_otu_table.otu_table_objects.append(list(genome_otu_table_collection))
            logging.info("Read in %i OTU sequences from reference genomes" % len(list(genomes_otu_table)))

            # read metagenome OTU samples and counts
            metagenome_sequence_to_sample_and_count = {}
            for otu in metagenome_otu_table:
                if otu.sequence not in metagenome_sequence_to_sample_and_count:
                    metagenome_sequence_to_sample_and_count[otu.sequence] = []
                metagenome_sequence_to_sample_and_count[otu.sequence].append(
                    SampleAndCount(otu.sample_name, otu.count))
                
            genome_otu_sequences = set([otu.sequence for otu in genomes_otu_table])
                    
            # make into a single OTU table collection
            total_collection = genomes_otu_table
            total_collection.add_otu_table_collection(metagenome_otu_table)
            
            # cluster and count as we go
            sample_name_to_appraisal = {}
            
            for cluster in Clusterer().cluster(total_collection, cluster_identity):
                # Determine if this cluster contains a genome representative
                cluster_represented = False
                
                for otu in cluster.otus:
                    if otu.sequence in genome_otu_sequences:
                        cluster_represented = True
                        break
                
                # Iterate through sequences, adding to num_found or num_not_found
                for otu in cluster.otus:
                    if otu.sequence in metagenome_sequence_to_sample_and_count: # possible that a genome was recovered but there is no OTU from the metagenome
                        for sample_and_count in metagenome_sequence_to_sample_and_count[otu.sequence]:
                            if sample_and_count.sample_name not in sample_name_to_appraisal:
                                res = AppraisalResult()
                                res.metagenome_sample_name = sample_and_count.sample_name
                                sample_name_to_appraisal[sample_and_count.sample_name] = res
                                
                            appraisal = sample_name_to_appraisal[sample_and_count.sample_name]
                            if cluster_represented:
                                appraisal.num_found += sample_and_count.count
                            else:
                                appraisal.num_not_found += sample_and_count.count
                                
            app = Appraisal()
            app.appraisal_results = sample_name_to_appraisal.values()
            return app

            
        
    def print_appraisal(self, appraisal,
                        output_io=sys.stdout,
                        accounted_for_otu_table=None,
                        unaccounted_for_otu_table=None):
        '''print the Appraisal object overview to STDOUT'''
        
        output_io.write("\t".join(['sample','num_found','num_not_found','percent_found'])+"\n")
        founds = []
        not_founds = []
        
        def print_sample(num_found, num_not_found, sample, mypercent=None):
            if mypercent:
                percent = mypercent
            elif num_found + num_not_found == 0:
                percent = 0.0
            else:
                percent = float(num_found)/(num_found+num_not_found) * 100
            output_io.write("\t".join([sample, str(num_found), str(num_not_found), "%2.1f" % percent])+"\n")
            
        def mean(l):
            return float(sum(l))/len(l) if len(l) > 0 else float('nan')
            
        for appraisal_result in appraisal.appraisal_results:
            print_sample(appraisal_result.num_found,
                         appraisal_result.num_not_found,
                         appraisal_result.metagenome_sample_name)
            founds.append(appraisal_result.num_found)
            not_founds.append(appraisal_result.num_not_found)

        print_sample(sum(founds), sum(not_founds), 'total')
        
        means = []
        for i, num_found in enumerate(founds):
            num_not_found = not_founds[i]
            means.append(float(num_found)/(num_found+num_not_found))
        print_sample("%2.1f" % mean(founds), "%2.1f" % mean(not_founds), 'average',
                     mypercent=mean(means)*100)
        
class SampleAndCount:
    def __init__(self, sample_name, count):
        self.sample_name = sample_name
        self.count = count
        
class AppraisalResult:
    num_found = 0
    num_not_found = 0
    metagenome_sample_name = None
    
class Appraisal:
    appraisal_results = None