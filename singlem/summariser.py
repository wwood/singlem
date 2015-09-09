import tempfile
import extern
from otu_table import OtuTable
from collections import OrderedDict
import logging

class Summariser:
    @staticmethod
    def summarise(**kwargs):
        '''Summarise an OTU table'''
        krona_output_prefix = kwargs.pop('krona_output_prefix', None)
        otu_table = kwargs.pop('otu_table', None)        
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        if krona_output_prefix:
            # prep the array
            gene_to_sample_to_taxonomy_to_count = {}
            for otu in OtuTable.each(open(otu_table)):
                if otu.marker not in gene_to_sample_to_taxonomy_to_count:
                    gene_to_sample_to_taxonomy_to_count[otu.marker] = OrderedDict()
                if otu.sample_name not in gene_to_sample_to_taxonomy_to_count[otu.marker]:
                    gene_to_sample_to_taxonomy_to_count[otu.marker][otu.sample_name] = OrderedDict()
                tax = '; '.join(otu.taxonomy_array() + [otu.sequence])
                if tax in gene_to_sample_to_taxonomy_to_count[otu.marker][otu.sample_name]:
                    raise Exception("Unexpected duplicated sequence/taxonomy found in OTU table")
                gene_to_sample_to_taxonomy_to_count[otu.marker][otu.sample_name][tax] = otu.coverage
                
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
                        f.write('\t'.join([str(coverage)]+taxonomy.split('; ')))
                        f.write('\n')
                    f.flush()
                    cmd += ' %s,%s' % (f.name, sample)
                extern.run(cmd)
                for f in sample_name_to_tempfile.values():
                    f.close()
            
        else:
            raise Exception()
        
        
        