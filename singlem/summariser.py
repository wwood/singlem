import tempfile
import extern
from collections import OrderedDict
import logging
from otu_table import OtuTable

class Summariser:
    @staticmethod
    def summarise(**kwargs):
        '''Summarise an OTU table'''
        krona_output_prefix = kwargs.pop('krona_output_prefix')
        table_collection = kwargs.pop('table_collection')        
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        # prep the array
        gene_to_sample_to_taxonomy_to_count = Summariser._collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection)
            
        # write the input krona files
        sample_name_to_tempfile = OrderedDict()
        for gene, sample_to_taxonomy_to_count in gene_to_sample_to_taxonomy_to_count.iteritems():
            krona_output_file = '%s.%s.krona.html' % (krona_output_prefix, gene)
            logging.info("Writing krona %s" % krona_output_file)
            cmd = 'ktImportText -o %s' % krona_output_file
            for sample, taxonomy_to_count in sample_to_taxonomy_to_count.iteritems():
                f = tempfile.NamedTemporaryFile(prefix='singlem_for_krona')
                sample_name_to_tempfile[sample] = f
                
                for taxonomy, coverage in taxonomy_to_count.iteritems():
                    tax_split = taxonomy.split('; ')
                    if tax_split[0] == 'Root' and len(tax_split) > 1: tax_split = tax_split[1:]
                    f.write('\t'.join([str(coverage)]+tax_split))
                    f.write('\n')
                f.flush()
                cmd += ' %s,%s' % (f.name, sample)
            extern.run(cmd)
            for f in sample_name_to_tempfile.values():
                f.close()
                
    @staticmethod
    def _collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection,
                                                                     add_sequence_to_taxonomy=True,
                                                                     use_coverage=True):
        gene_to_sample_to_taxonomy_to_count = {}
        for otu in table_collection:
            if otu.marker not in gene_to_sample_to_taxonomy_to_count:
                gene_to_sample_to_taxonomy_to_count[otu.marker] = OrderedDict()
            if otu.sample_name not in gene_to_sample_to_taxonomy_to_count[otu.marker]:
                gene_to_sample_to_taxonomy_to_count[otu.marker][otu.sample_name] = OrderedDict()
            if add_sequence_to_taxonomy:
                tax = '; '.join(otu.taxonomy_array() + [otu.sequence])
            else:
                tax = otu.taxonomy
            if add_sequence_to_taxonomy and tax in gene_to_sample_to_taxonomy_to_count[otu.marker][otu.sample_name]:
                raise Exception("Unexpected duplicated sequence/taxonomy found in OTU table")
            if use_coverage:
                record = otu.coverage
            else:
                record = otu.count
            gene_to_sample_to_taxonomy_to_count[otu.marker][otu.sample_name][tax] = record
        return gene_to_sample_to_taxonomy_to_count
        
    @staticmethod
    def write_unifrac_format_file(**kwargs):
        '''Summarise an OTU table as 3 column tab separated OTU table:
        sample, taxonomy, count
        
        When >1 OTU maps to the same taxonomy, collapse them into one OTU,
        adding their coverages
        '''
        unifrac_output_prefix = kwargs.pop('unifrac_output_prefix')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        # prep the array
        gene_to_sample_to_taxonomy_to_count = Summariser._collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection, add_sequence_to_taxonomy=False, use_coverage=False)
        
        for gene, sample_to_taxonomy_to_count in gene_to_sample_to_taxonomy_to_count.items():
            unifrac_output = '%s.%s.unifrac' % (unifrac_output_prefix, gene)
            logging.info("Writing %s" % unifrac_output)
            with open(unifrac_output, 'w') as f:
                for sample, taxonomy_to_count in sample_to_taxonomy_to_count.items():
                    for taxonomy, count in taxonomy_to_count.items():
                        if taxonomy == '': continue #ignore unclassified
                        f.write("\t".join([taxonomy, sample, str(count)])+"\n")
        logging.info("Finished")
        
    @staticmethod
    def write_otu_table(**kwargs):
        output_table_io = kwargs.pop('output_table_io')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        OtuTable.write_otus_to(table_collection, output_table_io)
            
        
                        
        
        
        
        