import logging
import sys

from otu_table import OtuTable
from otu_table_collection import OtuTableCollection
from sequence_searcher import SequenceSearcher

class Appraiser:
    def appraise(self, **kwargs):
        '''Given a collection of OTU tables derived from samples, and OTU
        table(s) corresponding to a collection of recovered genomes, how
        much of the community has been recovered in those genomes?

        Parameters
        ----------
        kwargs:
            sequence_identity: float for 'near enough', None when an exact match is required.
        
        Returns
        -------
        An Appraisal object containing appraisals for each metagenome
        '''
        genome_otu_table_collection = kwargs.pop('genome_otu_table_collection')
        metagenome_otu_table_collection = kwargs.pop('metagenome_otu_table_collection')
        sequence_identity = kwargs.pop('sequence_identity', None)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        logging.info("Read in %i markers from the different genomes" %\
                     len(genome_otu_table_collection))
        filtered_genome_otus = \
            list(genome_otu_table_collection.excluded_duplicate_distinct_genes())
        logging.info("After excluding duplicate markers that may indicate "
                     "contamination, found %i markers" % len(filtered_genome_otus))
        
        if sequence_identity is None:
            genome_otu_sequences = set()
            genome_names = set()
            for otu in filtered_genome_otus:
                genome_otu_sequences.add(otu.sequence)
                genome_names.add(otu.sample_name)
            logging.info("Read in %i unique sequences from the %i reference genomes" %\
                         (len(genome_otu_sequences), len(genome_names)))
            
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
                    appraisal.found_otus.append(otu)
                else:
                    appraisal.num_not_found += count
                    appraisal.not_found_otus.append(otu)
                    
            app = Appraisal()
            app.appraisal_results = sample_name_to_appraisal.values()
            return app
        
        else:
            sample_name_to_appraisal = {}
            seen_otus = set()
            genome_otu_table = OtuTable()
            genome_otu_table.add(filtered_genome_otus)
            filtered_collection = OtuTableCollection()
            filtered_collection.otu_table_objects = [genome_otu_table]
            for uc in SequenceSearcher().global_search(metagenome_otu_table_collection,
                                             filtered_collection,
                                             sequence_identity):
                q = uc.query
                key = str([q.sample_name, q.sequence])
                if key in seen_otus:
                    logging.warn("Double-saw an OTU..")
                    continue
                else:
                    seen_otus.add(key)
                if q.sample_name not in sample_name_to_appraisal:
                    res = AppraisalResult()
                    res.metagenome_sample_name = q.sample_name
                    sample_name_to_appraisal[q.sample_name] = res
                    
                appraisal = sample_name_to_appraisal[q.sample_name]
                if uc.target is None:
                    appraisal.num_not_found += q.count
                    appraisal.not_found_otus.append(q)
                else:
                    appraisal.num_found += q.count
                    appraisal.found_otus.append(q)
                    
            app = Appraisal()
            app.appraisal_results = sample_name_to_appraisal.values()
            return app

            
        
    def print_appraisal(self, appraisal,
                        output_io=sys.stdout,
                        accounted_for_otu_table_io=None,
                        unaccounted_for_otu_table_io=None):
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
        
        if accounted_for_otu_table_io:
            accounted_for_table = OtuTable()
        if unaccounted_for_otu_table_io:
            unaccounted_for_table = OtuTable()
            
        for appraisal_result in appraisal.appraisal_results:
            print_sample(appraisal_result.num_found,
                         appraisal_result.num_not_found,
                         appraisal_result.metagenome_sample_name)
            founds.append(appraisal_result.num_found)
            not_founds.append(appraisal_result.num_not_found)
            if accounted_for_otu_table_io:
                accounted_for_table.add(appraisal_result.found_otus)
            if accounted_for_otu_table_io:
                unaccounted_for_table.add(appraisal_result.not_found_otus)
            
        print_sample(sum(founds), sum(not_founds), 'total')
        
        means = []
        for i, num_found in enumerate(founds):
            num_not_found = not_founds[i]
            means.append(float(num_found)/(num_found+num_not_found))
        print_sample("%2.1f" % mean(founds), "%2.1f" % mean(not_founds), 'average',
                     mypercent=mean(means)*100)
        
        if accounted_for_otu_table_io:
            accounted_for_table.write_to(accounted_for_otu_table_io)
        if unaccounted_for_otu_table_io:
            unaccounted_for_table.write_to(unaccounted_for_otu_table_io)

        
class AppraisalResult:
    num_found = 0
    num_not_found = 0
    metagenome_sample_name = None
    
    def __init__(self):
        self.found_otus = []
        self.not_found_otus = []
    
class Appraisal:
    appraisal_results = None
