import logging

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
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
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
        
    def print_appraisal(self, appraisal):
        '''print the Appraisal object overview to STDOUT'''
        
        print "\t".join(['sample','num_found','num_not_found','percent_found'])
        total_found = 0
        total_not_found = 0
        
        def print_sample(num_found, num_not_found, sample):
            if num_found + num_not_found == 0:
                percent = 0.0
            else:
                percent = float(num_found)/(num_found+num_not_found) * 100
            print "\t".join([sample, num_found, num_not_found, percent])
            
        for sample, appraisal in sample_name_to_appraisal.items():
            print_sample(appraisal.num_found, appraisal.num_not_found, sample)
            total_found += appraisal.num_found
            total_not_found += appraisal.num_not_found
        print_sample(total_found, total_not_found, 'total')
        
        
class AppraisalResult:
    num_found = 0
    num_not_found = 0
    metagenome_sample_name = None
    
class Appraisal:
    appraisal_results = None