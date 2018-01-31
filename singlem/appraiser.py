import logging
import sys

from otu_table import OtuTable
from otu_table_collection import OtuTableCollection
from sequence_searcher import SequenceSearcher
from appraisal_result import Appraisal, AppraisalResult

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
        sequence_identity = kwargs.pop('sequence_identity', None)
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

        if sequence_identity is None:
            genome_otu_sequences = set()
            genome_names = set()
            if appraising_binning:
                for otu in filtered_genome_otus:
                    genome_otu_sequences.add(otu.sequence)
                    genome_names.add(otu.sample_name)
            if appraising_assembly:
                assembly_sequences = set()
                for otu in assembly_otu_table_collection:
                    assembly_sequences.add(otu.sequence)
            if appraising_binning:
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
                if appraising_binning and otu.sequence in genome_otu_sequences:
                    appraisal.num_binned += count
                    appraisal.binned_otus.append(otu)
                    # Probably this 'if' condition is not necessary, but just to check.
                    if appraising_assembly and otu.sequence in assembly_sequences:
                        appraisal.num_assembled += count
                        appraisal.assembled_otus.append(otu)
                elif appraising_assembly and otu.sequence in assembly_sequences:
                    appraisal.num_assembled += count
                    appraisal.assembled_otus.append(otu)
                else:
                    appraisal.num_not_found += count
                    appraisal.not_found_otus.append(otu)

            app = Appraisal()
            app.appraisal_results = sample_name_to_appraisal.values()
            return app

        else:
            if appraising_binning:
                sample_to_binned = self._appraise_inexactly(
                    metagenome_otu_table_collection,
                    filtered_genome_otus,
                    sequence_identity)
                sample_to_building_block = sample_to_binned
            if assembly_otu_table_collection:
                sample_to_assembled = self._appraise_inexactly(
                    metagenome_otu_table_collection,
                    assembly_otu_table_collection,
                    sequence_identity)
                sample_to_building_block = sample_to_assembled

            app = Appraisal()
            app.appraisal_results = []
            for sample in sample_to_building_block.keys():
                res = AppraisalResult()
                res.metagenome_sample_name = sample
                seen_otu_sequences = set()
                if appraising_binning:
                    binning_building_block = sample_to_binned[sample]
                    res.num_binned = binning_building_block.num_found
                    res.binned_otus = binning_building_block.found_otus
                    for otu in res.binned_otus:
                        seen_otu_sequences.add(otu.sequence)
                if appraising_assembly:
                    assembled_building_block = sample_to_assembled[sample]
                    res.num_assembled = assembled_building_block.num_found
                    res.assembled_otus = assembled_building_block.found_otus
                    for otu in res.assembled_otus:
                        seen_otu_sequences.add(otu.sequence)
                not_seen_otus = []
                for otu in metagenome_otu_table_collection:
                    if otu.sample_name == sample and otu.sequence not in seen_otu_sequences:
                        not_seen_otus.append(otu)
                res.not_found_otus = not_seen_otus
                for otu in not_seen_otus:
                    res.num_not_found += otu.count
                app.appraisal_results.append(res)
            return app


    def _appraise_inexactly(self, metagenome_otu_table_collection,
                            found_otu_collection,
                            sequence_identity):
        '''Given a metagenome sample collection and OTUs 'found' either by binning or
        assembly, return a AppraisalBuildingBlock representing the OTUs that
        have been found, using inexact matching.

        '''
        found_otu_table = OtuTable()
        found_otu_table.add(found_otu_collection)
        found_collection = OtuTableCollection()
        found_collection.otu_table_objects = [found_otu_table]

        sample_to_building_block = {}

        for uc in SequenceSearcher().global_search(metagenome_otu_table_collection,
                                         found_otu_collection,
                                         sequence_identity):
            q = uc.query
            if q.sample_name in sample_to_building_block:
                appraisal = sample_to_building_block[q.sample_name]
            else:
                appraisal = AppraisalBuildingBlock()
                sample_to_building_block[q.sample_name] = appraisal

            if uc.target is not None:
                appraisal.num_found += q.count
                appraisal.found_otus.append(q)

        return sample_to_building_block



    def print_appraisal(self, appraisal,
                        doing_binning,
                        output_io=sys.stdout,
                        doing_assembly=False,
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
            if binned_otu_table_io:
                binned_table.add(appraisal_result.binned_otus)
            if unbinned_otu_table_io:
                unbinned_table.add(appraisal_result.assembled_not_binned_otus())
            if assembled_otu_table_io:
                assembled_table.add(appraisal_result.assembled_otus)
            if unaccounted_for_otu_table_io:
                unaccounted_for_table.add(appraisal_result.not_found_otus)

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

        if binned_otu_table_io:
            binned_table.write_to(binned_otu_table_io)
        if unbinned_otu_table_io:
            unbinned_table.write_to(unbinned_otu_table_io)
        if assembled_otu_table_io:
            assembled_table.write_to(assembled_otu_table_io)
        if unaccounted_for_otu_table_io:
            unaccounted_for_table.write_to(unaccounted_for_otu_table_io)


class AppraisalBuildingBlock:
    '''Can represent binned OTUs or assembled OTUs'''
    def __init__(self):
        self.num_found = 0
        self.found_otus = []
