import itertools
import tempfile
import extern
from collections import OrderedDict
import logging
import pandas
import Bio
import pandas as pd
import polars as pl

from .otu_table import OtuTable
from .rarefier import Rarefier
from .ordered_set import OrderedSet
from .archive_otu_table import ArchiveOtuTable
from .taxonomy import QUERY_BASED_ASSIGNMENT_METHOD, DIAMOND_ASSIGNMENT_METHOD, NO_ASSIGNMENT_METHOD
from .condense import CondensedCommunityProfile

class Summariser:
    @staticmethod
    def write_otu_table_krona(**kwargs):
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
        sample_to_gene_to_taxonomy_to_count = {}
        all_sample_names = set()
        all_gene_names = set()
        for gene, sample_to_taxonomy_to_count in gene_to_sample_to_taxonomy_to_count.items():
            all_gene_names.add(gene)
            for sample, taxonomy_to_count in sample_to_taxonomy_to_count.items():
                all_sample_names.add(sample)
                if sample not in sample_to_gene_to_taxonomy_to_count:
                    sample_to_gene_to_taxonomy_to_count[sample] = {}
                sample_to_gene_to_taxonomy_to_count[sample][gene] = taxonomy_to_count
        is_more_than_one_sample = len(sample_to_gene_to_taxonomy_to_count) > 1
        for sample in sorted(all_sample_names):
            for gene in sorted(all_gene_names):
                if gene in sample_to_gene_to_taxonomy_to_count[sample]:
                    f = tempfile.NamedTemporaryFile(prefix='singlem_for_krona',mode='w')
                    sample_tempfiles.append(f)

                    taxonomy_to_count = sample_to_gene_to_taxonomy_to_count[sample][gene]
                    for taxonomy, coverage in taxonomy_to_count.items():
                        tax_split = taxonomy.split('; ')
                        if tax_split[0] == 'Root' and len(tax_split) > 1: tax_split = tax_split[1:]
                        f.write('\t'.join([str(coverage)]+tax_split))
                        f.write('\n')
                    f.flush()
                    if is_more_than_one_sample:
                        display_name = '%s: %s' % (sample, gene)
                    else:
                        display_name = gene
                    cmd += " %s,'%s'" % (f.name, display_name)

        extern.run(cmd)
        for f in sample_tempfiles:
            f.close()

    @staticmethod
    def _collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection,
                                                                     add_sequence_to_taxonomy=True,
                                                                     use_sequence_as_taxonomy=False,
                                                                     use_coverage=True):
        if add_sequence_to_taxonomy and use_sequence_as_taxonomy:
            raise Exception("Cannot specify both add_sequence_to_taxonomy and use_sequence_as_taxonomy")
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
            elif use_sequence_as_taxonomy:
                tax = otu.sequence
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
    def write_unifrac_by_otu_format_file(**kwargs):
        '''Summarise an OTU table as 3 column tab separated OTU table:
        sample, OTU sequence, count

        When >1 OTU has the same sequence, collapse them into one OTU,
        adding their counts.
        '''
        unifrac_output_prefix = kwargs.pop('unifrac_output_prefix')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # prep the array
        gene_to_sample_to_taxonomy_to_count = Summariser._collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection, add_sequence_to_taxonomy=False, use_sequence_as_taxonomy=True, use_coverage=False)

        for gene, sample_to_taxonomy_to_count in gene_to_sample_to_taxonomy_to_count.items():
            unifrac_output = '%s.%s.unifrac' % (unifrac_output_prefix, gene)
            logging.info("Writing %s" % unifrac_output)
            with open(unifrac_output, 'w') as f:
                for sample, taxonomy_to_count in sample_to_taxonomy_to_count.items():
                    for taxonomy, count in taxonomy_to_count.items():
                        if taxonomy == '': continue #ignore unclassified
                        f.write("\t".join([taxonomy, sample, str(count)])+"\n")

    @staticmethod
    def write_unifrac_by_taxonomy_format_file(**kwargs):
        '''Summarise an OTU table as 3 column tab separated OTU table:
        sample, taxonomy, count

        When >1 OTU has the same taxonomy, collapse them into one OTU,
        adding their counts.
        '''
        unifrac_output_prefix = kwargs.pop('unifrac_output_prefix')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # prep the array
        gene_to_sample_to_taxonomy_to_count = Summariser._collapse_otu_table_into_gene_to_sample_to_taxonomy_to_count(table_collection, add_sequence_to_taxonomy=False, use_sequence_as_taxonomy=False, use_coverage=False)

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
    def write_wide_format_otu_table(**kwargs):
        output_table_io = kwargs.pop('output_table_io')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if hasattr(output_table_io, 'name'):
            logging.info("Writing %s" % output_table_io.name)
        else:
            logging.info("Writing an OTU table")

        # Collect a hash of sequence to sample to num_seqs
        gene_to_seq_to_sample_to_count = OrderedDict()
        sequence_to_taxonomy = {}
        samples = OrderedSet()
        for otu in table_collection:
            if otu.marker not in gene_to_seq_to_sample_to_count:
                gene_to_seq_to_sample_to_count[otu.marker] = {}
            if otu.sequence not in gene_to_seq_to_sample_to_count[otu.marker]:
                gene_to_seq_to_sample_to_count[otu.marker][otu.sequence] = {}
            if otu.sample_name in gene_to_seq_to_sample_to_count[otu.marker][otu.sequence]:
                raise Exception("Unexpectedly found 2 of the same sequences for the same sample and marker")
            gene_to_seq_to_sample_to_count[otu.marker][otu.sequence][otu.sample_name] = otu.count
            samples.add(otu.sample_name)
            # This isn't perfect, because the same sequence might have
            # different taxonomies in different samples. But taxonomy might
            # be of regular form, or as a diamond example etc, so eh.
            sequence_to_taxonomy[otu.sequence] = otu.taxonomy

        output_table_io.write("\t".join(itertools.chain( # header
            ['marker','sequence'],
            samples,
            ['taxonomy\n'])))
        for gene, seq_to_sample_to_count in gene_to_seq_to_sample_to_count.items():
            for seq, sample_to_count in seq_to_sample_to_count.items():
                row = [gene, seq]
                for sample in samples:
                    try:
                        row.append(str(sample_to_count[sample]))
                    except KeyError:
                        row.append('0')
                row.append(sequence_to_taxonomy[seq])
                output_table_io.write("\t".join(row)+"\n")

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

    @staticmethod
    def write_translated_otu_table(**kwargs):
        output_table_io = kwargs.pop('output_table_io')
        table_collection = kwargs.pop('table_collection')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        printed_header = False
        def print_chunk(printed_header, seq_to_otus, output_table_io):
            if not printed_header:
                eg_otu = list(seq_to_otus.values())[0][0]
                output_table_io.write('\t'.join(
                    eg_otu.fields
                )+'\n')
            for (translated, to_collapse) in seq_to_otus.items():
                total_num_hits = 0
                total_coverage = 0
                max_hits = 0
                max_hits_taxonomy = 'programming error'
                for otu in to_collapse:
                    total_num_hits += otu.count
                    total_coverage += otu.coverage
                    if otu.count > max_hits:
                        max_hits_taxonomy = otu.taxonomy
                output_table_io.write('\t'.join([
                    otu.marker,
                    otu.sample_name,
                    translated,
                    str(total_num_hits),
                    str(total_coverage),
                    max_hits_taxonomy
                ])+'\n')

        original_count = 0
        collapsed_count = 0
        seq_to_otus = {}
        last_sample_and_marker = None
        for otu in table_collection:
            original_count += 1

            if last_sample_and_marker is None:
                last_sample_and_marker = [otu.sample_name, otu.marker]
            elif last_sample_and_marker != [otu.sample_name, otu.marker]:
                collapsed_count += len(seq_to_otus)
                print_chunk(printed_header, seq_to_otus, output_table_io)
                if not printed_header:
                    printed_header = True
                seq_to_otus = {}

            seq = str(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(otu.sequence)).translate().seq)
            if seq in seq_to_otus:
                seq_to_otus[seq].append(otu)
            else:
                seq_to_otus[seq] = [otu]

        collapsed_count += len(seq_to_otus)
        print_chunk(printed_header, seq_to_otus, output_table_io)
        logging.info("Printed {} collapsed OTUs from {} original OTUs".format(
            collapsed_count, original_count
        ))

    @staticmethod
    # args.collapse_paired_with_unpaired:
    # Summariser.write_collapsed_paired_with_unpaired_otu_table(
    #     archive_otu_tables = args.input_archive_otu_tables,
    #     output_table_io = open(args.collapse_paired_with_unpaired,'w'))
    def write_collapsed_paired_with_unpaired_otu_table(**kwargs):
        archive_otu_tables = kwargs.pop('archive_otu_tables')
        output_table_io = kwargs.pop('output_table_io')
        set_sample_name = kwargs.pop('set_sample_name', None)  # For merging OTU tables
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # Read all OTU tables
        df = None
        for a in archive_otu_tables:
            with open(a) as f:
                logging.debug("Reading archive table {} into RAM ..".format(a))
                ar = ArchiveOtuTable.read(f)
            if df is None:
                version = ar.version
                fields = ar.fields
                alignment_hmm_sha256s = ar.alignment_hmm_sha256s
                singlem_package_sha256s = ar.singlem_package_sha256s
                df = pandas.DataFrame(ar.data)
                df.columns = fields
            else:
                if version != ar.version:
                    raise Exception("Version mismatch between archives")
                elif fields != ar.fields:
                    raise Exception("Fields mismatch between archives")
                elif alignment_hmm_sha256s != ar.alignment_hmm_sha256s:
                    raise Exception("Alignment HMM SHA256 mismatch between archives")
                elif singlem_package_sha256s != ar.singlem_package_sha256s:
                    raise Exception("Singlem package SHA256 mismatch between archives")
                df2 = pandas.DataFrame(ar.data)
                df2.columns = fields
                df = pd.concat([df, df2], ignore_index=True)

        # Remove suffixes
        if set_sample_name is None:
            def remove_suffix(s):
                if s.endswith('_1'):
                    return s[:-2]
                else:
                    return s
            df['sample'] = df['sample'].apply(remove_suffix)

        # Ensure that there is now only exactly 1 sample name
        if set_sample_name is None and len(df['sample'].unique()) != 1:
            raise Exception("Multiple sample names found: {}".format(', '.join(df['sample'].unique())))
        if len(df['taxonomy_by_known?'].unique()) != 1:
            raise Exception("Multiple taxonomy_by_known found: {}".format(', '.join(df['taxonomy_by_known'].unique())))

        def combine_rows(grouped1):
            grouped = grouped1.reset_index()
            max_row = grouped['num_hits'].idxmax()
            if set_sample_name:
                sample = set_sample_name
            else:
                sample = grouped.iloc[0]['sample']
            tax_assignment_method = grouped.iloc[0]['taxonomy_assignment_method']
            if tax_assignment_method == QUERY_BASED_ASSIGNMENT_METHOD:
                equal_best_hit_taxonomies = grouped.iloc[0]['equal_best_hit_taxonomies']
            elif tax_assignment_method == DIAMOND_ASSIGNMENT_METHOD:
                equal_best_hit_taxonomies = list(itertools.chain(*grouped['equal_best_hit_taxonomies']))
            elif tax_assignment_method == None or tax_assignment_method == NO_ASSIGNMENT_METHOD:
                equal_best_hit_taxonomies = None
            else:
                raise Exception("Unexpected tax assignment method: {}".format(tax_assignment_method))
            return pd.DataFrame({
                'gene':[grouped.iloc[0]['gene']],
                'sample':[sample],
                'sequence':[grouped.iloc[0]['sequence']],
                'num_hits':[sum(grouped['num_hits']),],
                'coverage':[sum(grouped['coverage']),],
                'taxonomy':[grouped.iloc[max_row]['taxonomy']],
                'read_names':[list(itertools.chain(*grouped['read_names']))],
                'nucleotides_aligned':[list(itertools.chain(*grouped['nucleotides_aligned']))],
                'taxonomy_by_known?':[grouped.iloc[0]['taxonomy_by_known?']],
                'read_unaligned_sequences':[list(itertools.chain(*grouped['read_unaligned_sequences']))],
                'equal_best_hit_taxonomies':[equal_best_hit_taxonomies],
                'taxonomy_assignment_method':[tax_assignment_method],
            })
        transformed = df.groupby(['sequence','gene'], as_index=False).apply(combine_rows)[ArchiveOtuTable.FIELDS]
        logging.info("Collapsed {} total OTUs into {} output OTUs".format(len(df), len(transformed)))

        logging.debug("Writing output table ..")
        ar.data = transformed.values.tolist()
        ar.write_to(output_table_io)
        # json.dump({"version": ar.version,
        #     "alignment_hmm_sha256s": ar.alignment_hmm_sha256s,
        #     "singlem_package_sha256s": ar.singlem_package_sha256s,
        #     'fields': ar.fields,
        #     "otus": transformed.to_json(orient='values')},
        #     output_table_io)
        logging.info("Finished writing collapsed output table")

    @staticmethod
    def dump_raw_sequences_from_archive_otu_table(**kwargs):
        input_otus = kwargs.pop('table_collection')
        output_table_io = kwargs.pop('output_table_io')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        num_samples = 0
        num_otus = 0
        num_reads = 0
        seq_name_counters = {}
        for sample, otus in input_otus.each_sample_otus(generate_archive_otu_table=True):
            num_samples += 1
            for otu in otus:
                num_otus += 1
                for read_name, read_seq in zip(otu.read_names(), otu.read_unaligned_sequences()):
                    num_reads += 1
                    if read_name not in seq_name_counters:
                        seq_name_counters[read_name] = 0
                    output_table_io.write(">{}~{}\n{}\n".format(read_name, seq_name_counters[read_name], read_seq))
                    seq_name_counters[read_name] += 1
        logging.info("Wrote {} reads from {} OTUs in {} samples".format(num_reads, num_otus, num_samples))

    @staticmethod
    def write_species_by_site_table(**kwargs):
        input_taxonomic_profiles = kwargs.pop('input_taxonomic_profiles')
        output_species_by_site_relative_abundance_table = kwargs.pop('output_species_by_site_relative_abundance_table')
        level = kwargs.pop('level')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        logging.info("Writing site by species table")
        levels = ['root','domain','phylum','class','order','family','genus','species']
        if level not in levels:
            raise Exception("Unexpected level: %s" % level)
        level_index = levels.index(level)

        # Read the taxonomic profile
        all_profiles = []
        for profile_file in input_taxonomic_profiles:
            with open(profile_file) as f:
                for profile in CondensedCommunityProfile.each_sample_wise(f):
                    unassigned = 0.
                    name_to_coverage = {}
                    for node in profile.breadth_first_iter():
                        node_level = node.calculate_level()
                        if node_level < level_index:
                            unassigned += node.coverage
                        elif node_level == level_index:
                            name_to_coverage['; '.join(node.get_taxonomy())] = node.get_full_coverage()
                        else:
                            # node_level > level_index
                            # This coveage already included 
                            pass
                    # TODO: Sort keys? Later down I guess
                    keys = list(name_to_coverage.keys()) # not 100% confident that iterations will be in same order
                    current_sample = pl.DataFrame({
                        'taxonomy': keys,
                        'coverage': [name_to_coverage[k] for k in keys]
                    })
                    current_sample.columns = ['taxonomy','coverage']
                    current_sample = pl.concat([
                        pl.DataFrame([{'taxonomy':'unassigned','coverage':unassigned}]),
                        current_sample
                    ])
                    current_sample = current_sample.with_columns(pl.lit(profile.sample).alias('sample'))
                    # Make into a relabund
                    total_coverage = current_sample['coverage'].sum()
                    current_sample = current_sample.with_columns((pl.col('coverage') / total_coverage * 100).round(2).alias('relative_abundance'))
                    all_profiles.append(current_sample)

        # Concatenate all profiles
        all_profiles = pl.concat(all_profiles)
        all_profiles.pivot(index=['taxonomy'], columns=['sample'], values=['relative_abundance']).fill_null(0).write_csv(output_species_by_site_relative_abundance_table, separator='\t')
        logging.info("Wrote site by species table for {} samples".format(len(all_profiles['sample'].unique())))

    @staticmethod
    def write_taxonomic_level_coverage_table(**kwargs):
        input_taxonomic_profiles = kwargs.pop('input_taxonomic_profiles')
        output_taxonomic_level_coverage_table = kwargs.pop('output_taxonomic_level_coverage_table')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        logging.info("Writing taxonomic level coverage table")
        # Read the taxonomic profile
        all_profiles = []
        for profile_file in input_taxonomic_profiles:
            with open(profile_file) as f:
                for profile in CondensedCommunityProfile.each_sample_wise(f):
                    name_to_coverage = {}
                    for node in profile.breadth_first_iter():
                        node_level = node.calculate_level()
                        if node_level == 0:
                            continue
                        if node_level not in name_to_coverage:
                            name_to_coverage[node_level] = 0.
                        name_to_coverage[node_level] += node.coverage
                    all_profiles.append(pl.DataFrame({
                        'level': list(name_to_coverage.keys()),
                        'coverage': list(name_to_coverage.values())
                    }).with_columns(pl.lit(profile.sample).alias('sample')).with_columns(
                        ((pl.col('coverage') / pl.col('coverage').sum()).alias('relative_abundance') * 100).round(2),
                    ))
        # Concatenate all profiles
        all_profiles = pl.concat(all_profiles)

        if len(all_profiles.select(pl.col('level')).groupby('level').count()) in [7, 8]:
            # If there's 7 or 8 (including 0) levels, then assume that this is a regular taxonomy going on.
            levels = ['root','domain','phylum','class','order','family','genus','species']
            all_profiles = all_profiles.with_columns(pl.lit([levels[l] for l in all_profiles['level']]).alias('level'))

        all_profiles = all_profiles.select([
            'sample',
            'level',
            pl.col('coverage').round(2),
            pl.col('relative_abundance').alias('relative abundance (%)'),
        ])

        all_profiles.write_csv(output_taxonomic_level_coverage_table, separator='\t')

    def write_taxonomic_profile(input_taxonomic_profiles, output_taxonomic_profile_io):
        '''Write a taxonomic profile to a file'''
        logging.info("Writing taxonomic profile")
        seen_samples = set()
        CondensedCommunityProfile.write_header_to(output_taxonomic_profile_io)
        for profile_file in input_taxonomic_profiles:
            with open(profile_file) as f:
                for profile in CondensedCommunityProfile.each_sample_wise(f):
                    if profile.sample in seen_samples:
                        raise Exception("Duplicate sample name detected: %s" % profile.sample)
                    profile.write_data_to(output_taxonomic_profile_io)
                    seen_samples.add(profile.sample)
        logging.info("Wrote taxonomic profile")
        