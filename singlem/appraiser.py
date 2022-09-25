import logging
import os
import sys
import tempfile
import numpy

from .otu_table import OtuTable
from .otu_table_collection import OtuTableCollection
from .sequence_searcher import SequenceSearcher
from .appraisal_result import Appraisal, AppraisalResult
from .querier import Querier
from .sequence_database import SequenceDatabase

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
        assembly_otu_table_collection = kwargs.pop('assembly_otu_table_collection', None)
        packages = kwargs.pop('packages')
        sequence_identity = kwargs.pop('sequence_identity', None)
        output_found_in = kwargs.pop('output_found_in', False)
        window_size = kwargs.pop('window_size')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        appraising_binning = genome_otu_table_collection is not None
        if appraising_binning:
            logging.info("Read in %i markers from the different genomes" %\
                         len(genome_otu_table_collection))
            filtered_genome_otus = \
                list(genome_otu_table_collection.excluded_duplicate_distinct_genes())
            logging.info("After excluding duplicate markers that may indicate "
                         "contamination, found %i markers" % len(filtered_genome_otus))
            if len(filtered_genome_otus) == 0:
                logging.warning("No markers suitable for appraisal found. This may be because all found markers "
                            "have been excluded. To run 'appraise', the source genome name should be "
                            "reflected in the 'sample' column of the OTU table, i.e. singlem 'pipe' should "
                            "be run on several FASTA files, not a concatenation of genomes.")

        appraising_assembly = assembly_otu_table_collection is not None

        if appraising_binning:
            sample_to_binned = self._appraise_inexactly(
                metagenome_otu_table_collection,
                filtered_genome_otus,
                sequence_identity,
                output_found_in,
                packages,
                window_size)
            sample_to_building_block = sample_to_binned
        if assembly_otu_table_collection:
            sample_to_assembled = self._appraise_inexactly(
                metagenome_otu_table_collection,
                assembly_otu_table_collection,
                sequence_identity,
                output_found_in,
                packages,
                window_size)
            sample_to_building_block = sample_to_assembled

        app = Appraisal()
        app.appraisal_results = []
        for sample in list(sample_to_building_block.keys()):
            res = AppraisalResult()
            res.metagenome_sample_name = sample
            seen_otu_sequences = set()
            if appraising_binning:
                binning_building_block = sample_to_binned[sample]
                res.num_binned = binning_building_block.est_num_found()
                res.binned_otus = binning_building_block.found_otus
                for otu in res.binned_otus:
                    seen_otu_sequences.add(otu.sequence)
            if appraising_assembly:
                assembled_building_block = sample_to_assembled[sample]
                res.num_assembled = assembled_building_block.est_num_found()
                res.assembled_otus = assembled_building_block.found_otus
                for otu in res.assembled_otus:
                    seen_otu_sequences.add(otu.sequence)
            not_seen_otus = []
            for otu in metagenome_otu_table_collection:
                if otu.sample_name == sample and otu.sequence not in seen_otu_sequences:
                    if output_found_in:
                        otu.add_found_data('')
                    not_seen_otus.append(otu)
            res.not_found_otus = not_seen_otus
            not_found_block = AppraisalBuildingBlock(packages)
            for otu in not_seen_otus:
                not_found_block.add_otu(otu)
            res.num_not_found = not_found_block.est_num_found()
            app.appraisal_results.append(res)
        return app


    def _appraise_inexactly(self, metagenome_otu_table_collection,
                            found_otu_collection,
                            sequence_identity,
                            output_found_in,
                            packages,
                            window_size):
        '''Given a metagenome sample collection and OTUs 'found' either by binning or
        assembly, return a AppraisalBuildingBlock representing the OTUs that
        have been found, using inexact matching.

        '''
        if sequence_identity:
            max_divergence = window_size * (1 - sequence_identity)
        else:
            max_divergence = 0

        tmp = tempfile.TemporaryDirectory()
        sdb_path = os.path.join(tmp.name, "tmp.sdb")
        sequence_database = SequenceDatabase()
        sequence_database.create_from_otu_table(sdb_path, found_otu_collection, sequence_database_methods = ['naive'])
        sdb_tmp = sequence_database.acquire(sdb_path)

        found_genes = [table.marker for table in found_otu_collection]
        metagenome_table = OtuTable()
        for otu in metagenome_otu_table_collection:
            if otu.marker in found_genes:
                metagenome_table.add([otu])

        metagenome_collection = OtuTableCollection()
        metagenome_collection.add_otu_table_object(metagenome_table)
        metagenome_collection.sort_otu_tables_by_marker()

        querier = Querier()
        queries = querier.query_with_queries(metagenome_collection, sdb_tmp, max_divergence, 'naive', SequenceDatabase.NUCLEOTIDE_TYPE, 1, None, True, None)

        sample_to_building_block = {}
        for hit in queries:
            # hit has (query, subject, divergence)
            # subject has .taxonomy
            q = hit.query
            if q.sample_name in sample_to_building_block:
                appraisal = sample_to_building_block[q.sample_name]
            else:
                appraisal = AppraisalBuildingBlock(packages)
                sample_to_building_block[q.sample_name] = appraisal

            if output_found_in:
                q.add_found_data(hit.subject.sample_name)

            if q not in appraisal.found_otus:
                appraisal.add_otu(q)

        for otu in metagenome_otu_table_collection:
            if otu.sample_name not in sample_to_building_block:
                sample_to_building_block[otu.sample_name] = AppraisalBuildingBlock(packages)

        return sample_to_building_block



    def print_appraisal(self, appraisal,
                        doing_binning,
                        output_io=sys.stdout,
                        doing_assembly=False,
                        output_found_in=False,
                        binned_otu_table_io=None,
                        unbinned_otu_table_io=None,
                        assembled_otu_table_io=None,
                        unaccounted_for_otu_table_io=None):
        '''print the Appraisal object overview to STDOUT'''

        headers = ['sample']
        if doing_binning: headers.append('num_binned')
        if doing_assembly: headers.append('num_assembled')
        headers.append('num_not_found')
        if doing_binning: headers.append('percent_binned')
        if doing_assembly: headers.append('percent_assembled')
        output_io.write("\t".join(headers)+"\n")

        binned = []
        assembled = []
        assembled_not_binned = []
        not_founds = []

        def print_sample(num_binned, num_assembled, num_assembled_not_binned, num_not_found, sample,
                         mypercent_binned=None, mypercent_assembled=None):
            if mypercent_binned is not None or mypercent_assembled is not None:
                if doing_binning:
                    percent_binned = mypercent_binned
                if doing_assembly:
                    percent_assembled = mypercent_assembled
            else:
                total = num_not_found
                if doing_binning: total += num_binned
                if doing_assembly: total += num_assembled_not_binned
                if total == 0:
                    if doing_binning: percent_binned = 0.0
                    if doing_assembly: percent_assembled = 0.0
                else:
                    if doing_binning:
                        percent_binned = float(num_binned)/total * 100
                    if doing_assembly:
                        percent_assembled = float(num_assembled)/total * 100
            to_write = [sample]
            if doing_binning: to_write.append(str(num_binned))
            if doing_assembly: to_write.append(str(num_assembled))
            to_write.append(str(num_not_found))
            if doing_binning:
                to_write.append("%2.1f" % percent_binned)
            if doing_assembly:
                to_write.append("%2.1f" % percent_assembled)
            output_io.write("\t".join(to_write)+"\n")

        def mean(l):
            return float(sum(l))/len(l) if len(l) > 0 else float('nan')

        if binned_otu_table_io:
            binned_table = OtuTable()
        if unbinned_otu_table_io:
            unbinned_table = OtuTable()
        if assembled_otu_table_io:
            assembled_table = OtuTable()
        if unaccounted_for_otu_table_io:
            unaccounted_for_table = OtuTable()

        for appraisal_result in appraisal.appraisal_results:
            if doing_assembly:
                num_assembled_not_binned = appraisal_result.num_assembled_not_binned()
            print_sample(appraisal_result.num_binned if doing_binning else None,
                         appraisal_result.num_assembled if doing_assembly else None,
                         num_assembled_not_binned if doing_assembly else None,
                         appraisal_result.num_not_found,
                         appraisal_result.metagenome_sample_name)
            if doing_binning:
                binned.append(appraisal_result.num_binned)
            if doing_assembly:
                assembled.append(appraisal_result.num_assembled)
                assembled_not_binned.append(num_assembled_not_binned)
            not_founds.append(appraisal_result.num_not_found)

            if not output_found_in:
                if binned_otu_table_io:
                    binned_table.add(appraisal_result.binned_otus)
                if unbinned_otu_table_io:
                    unbinned_table.add(appraisal_result.assembled_not_binned_otus())
                if assembled_otu_table_io:
                    assembled_table.add(appraisal_result.assembled_otus)
                if unaccounted_for_otu_table_io:
                    unaccounted_for_table.add(appraisal_result.not_found_otus)
            else:
                if binned_otu_table_io:
                    binned_table.add_with_extras(appraisal_result.binned_otus, ['found_in'])
                if unbinned_otu_table_io:
                    unbinned_table.add_with_extras(appraisal_result.assembled_not_binned_otus(), ['found_in'])
                if assembled_otu_table_io:
                    assembled_table.add_with_extras(appraisal_result.assembled_otus, ['found_in'])
                if unaccounted_for_otu_table_io:
                    unaccounted_for_table.add_with_extras(appraisal_result.not_found_otus, ['found_in'])

        print_sample(sum(binned) if doing_binning else None,
                     sum(assembled) if doing_assembly else None,
                     sum(assembled_not_binned) if doing_assembly else None,
                     sum(not_founds),
                     'total')

        binned_means = []
        assembled_means = []
        if doing_binning:
            to_enumerate = binned
        else:
            to_enumerate = assembled
        for i, _ in enumerate(to_enumerate):
            num_binned = binned[i] if doing_binning else 0
            num_assembled = assembled[i] if doing_assembly else 0
            num_assembled_not_binned = assembled_not_binned[i] if doing_assembly else 0
            num_not_found = not_founds[i]
            total = num_assembled_not_binned+num_not_found
            if doing_binning:
                total += num_binned
                binned_means.append(float(num_binned)/total)
            if doing_assembly:
                assembled_means.append(float(num_assembled)/total)
        print_sample("%2.1f" % mean(binned) if doing_binning else None,
                     "%2.1f" % mean(assembled) if doing_assembly else None,
                     None,
                     "%2.1f" % mean(not_founds),
                     'average',
                     mypercent_binned=mean(binned_means)*100 if doing_binning else None,
                     mypercent_assembled=(mean(assembled_means)*100 if doing_assembly else None))

        if not output_found_in:
            if binned_otu_table_io:
                binned_table.write_to(binned_otu_table_io)
            if unbinned_otu_table_io:
                unbinned_table.write_to(unbinned_otu_table_io)
            if assembled_otu_table_io:
                assembled_table.write_to(assembled_otu_table_io)
            if unaccounted_for_otu_table_io:
                unaccounted_for_table.write_to(unaccounted_for_otu_table_io)
        else:
            if binned_otu_table_io:
                binned_table.write_to(binned_otu_table_io, fields_to_print=binned_table.fields)
            if unbinned_otu_table_io:
                unbinned_table.write_to(unbinned_otu_table_io, fields_to_print=binned_table.fields)
            if assembled_otu_table_io:
                assembled_table.write_to(assembled_otu_table_io, fields_to_print=binned_table.fields)
            if unaccounted_for_otu_table_io:
                unaccounted_for_table.write_to(unaccounted_for_otu_table_io, fields_to_print=binned_table.fields)


class AppraisalBuildingBlock:
    '''Can represent binned OTUs or assembled OTUs'''
    DOMAINS = ['d__Archaea', 'd__Bacteria']

    def __init__(self, packages):
        self.num_found = {}
        for package in packages:
            domains = ['d__' + d for d in package.target_domains()]
            self.num_found[package.graftm_package_basename()] = {domain: 0 for domain in domains if domain in self.DOMAINS}

        self.found_otus = []

    def add_otu(self, otu):
        if otu.marker in self.num_found:
            try:
                domain = otu.taxonomy.split("; ")[1]
                self.num_found[otu.marker][domain] += otu.count
            except (IndexError, KeyError):
                # If the taxonomy is not defined, or the domain is not in target_domain
                pass
        else:
            raise Exception("Unexpected marker found: %s" % otu.marker)
        self.found_otus.append(otu)

    def est_num_found(self):
        out = {}
        for domain in self.DOMAINS:
            out[domain] = round(self._trimmean([n[domain] for n in self.num_found.values() if domain in n], 10))
        return out

    def _trimmean(self, arr, percent):
        n = len(arr)
        k = int(round(n*(float(percent)/100)/2))
        return numpy.mean(arr[k+1:n-k])
