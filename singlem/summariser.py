import tempfile
import extern
from collections import OrderedDict
import logging
from otu_table import OtuTable
from rarefier import Rarefier

class Summariser:
    @staticmethod
    def summarise(**kwargs):
        '''Summarise an OTU table'''
        krona_output_file = kwargs.pop('krona_output')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # prep the array
        gene_to_sample_to_taxonomy_to_count = Summariser._collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection)

        # write the output krona files
        sample_name_to_tempfile = OrderedDict()
        logging.info("Writing krona %s" % krona_output_file)
        cmd = 'ktImportText -o %s' % krona_output_file
        sample_tempfiles = []
        for gene, sample_to_taxonomy_to_count in gene_to_sample_to_taxonomy_to_count.iteritems():
            for sample, taxonomy_to_count in sample_to_taxonomy_to_count.iteritems():
                f = tempfile.NamedTemporaryFile(prefix='singlem_for_krona')
                sample_tempfiles.append(f)

                for taxonomy, coverage in taxonomy_to_count.iteritems():
                    tax_split = taxonomy.split('; ')
                    if tax_split[0] == 'Root' and len(tax_split) > 1: tax_split = tax_split[1:]
                    f.write('\t'.join([str(coverage)]+tax_split))
                    f.write('\n')
                f.flush()
                if len(sample_to_taxonomy_to_count) == 1:
                    display_name = gene
                else:
                    display_name = '%s: %s' % (sample, gene)
                cmd += " %s,'%s'" % (f.name, display_name)
        extern.run(cmd)
        for f in sample_tempfiles:
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
                if otu.taxonomy_array():
                    tax = '; '.join(otu.taxonomy_array() + [otu.sequence])
                else:
                    tax = '; '.join([otu.sequence])
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

    @staticmethod
    def write_otu_table(**kwargs):
        output_table_io = kwargs.pop('output_table_io')
        table_collection = kwargs.pop('table_collection')
        output_extras = kwargs.pop('output_extras')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if hasattr(output_table_io, 'name'):
            logging.info("Writing %s" % output_table_io.name)
        else:
            logging.info("Writing an OTU table")

        if output_extras:
            OtuTable.write_otus_to(table_collection, output_table_io,
                                   fields_to_print=table_collection.example_field_names())
        else:
            OtuTable.write_otus_to(table_collection, output_table_io)

    @staticmethod
    def write_clustered_otu_table(**kwargs):
        output_table_io = kwargs.pop('output_table_io')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        logging.info("Writing clustered OTU table")
        output_table_io.write(
            "\t".join(
                OtuTable.DEFAULT_OUTPUT_FIELDS+
                ['representative',
                 'total_num_reads',
                 'total_coverage',
                 'num_sub_otus',
                 'max_sub_otu_abundance'])
            +"\n")

        for d in table_collection:
            for otu in d.otus:
                output_table_io.write("\t".join(
                    [OtuTable._to_printable(cell) for cell in [
                        otu.marker,
                        otu.sample_name,
                        otu.sequence,
                        otu.count,
                        otu.coverage,
                        otu.taxonomy,
                        d.sequence,
                        d.count,
                        d.coverage,
                        len(d.otus),
                        max([otu.count for otu in d.otus])
                    ]])+"\n")

    @staticmethod
    def write_rarefied_otu_table(**kwargs):
        output_table_io = kwargs.pop('output_table_io')
        table_collection = kwargs.pop('table_collection')
        number_to_choose = kwargs.pop('number_to_choose', None)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if number_to_choose is None:
            counts = {}
            for otu in table_collection:
                key = "%s_singlem_RAND8_%s" % (otu.sample_name, otu.marker)
                try:
                    counts[key] += otu.count
                except KeyError:
                    counts[key] = otu.count
            number_to_choose = min(counts.values())
            logging.info("Minimum number of sequences detected is %i, rarefying all sample/gene combinations to this level" % number_to_choose)

        logging.info("Rarefying OTU table to max %i sequences per sample/gene combination and writing to %s" % (number_to_choose, output_table_io.name))
        OtuTable.write_otus_to(Rarefier().rarefy(table_collection, number_to_choose),
                               output_table_io)
