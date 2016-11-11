import logging

from archive_otu_table import ArchiveOtuTable
from otu_table import OtuTable
from taxonomy import Taxonomy

class OtuTableCollection:
    def __init__(self):
        self.otu_table_objects = []
        self.archive_table_objects = []

        # None or an of taxonomy to iterate over
        self.target_taxonomy = None

    def add_otu_table(self, input_otu_table_io):
        '''Add a regular style OTU table to the collection.

        Parameters
        ----------
        input_otu_table_ios: list of IO
            entries are open streams of OTU table data

        Returns
        -------
        None
        '''
        self.otu_table_objects.append(OtuTable.read(input_otu_table_io))

    def add_archive_otu_table(self, input_archive_table_io):
        self.archive_table_objects.append(ArchiveOtuTable.read(input_archive_table_io))

    def add_otu_table_collection(self, otu_table_collection):
        '''Append an OtuTableCollection to this collection.
        Only the tables are added, the target_taxonomy is ignored'''
        for otu_table in otu_table_collection.otu_table_objects:
            self.otu_table_objects.append(otu_table)
        for archive in otu_table_collection.archive_table_objects:
            self.archive_table_objects.append(archive)

    def set_target_taxonomy_by_string(self, taxonomy_string):
        '''Set the target_taxonomy instance variable by a string, which
        gets parsed into the requisite array form and stored in the instance
        variable'''
        self.target_taxonomy = Taxonomy.split_taxonomy(taxonomy_string)

    def example_field_names(self):
        '''Return the field names of the first OTU table'''
        for table_types in (self.otu_table_objects, self.archive_table_objects):
            for table in table_types:
                return table.fields
        raise Exception("Attempt to get fields from empty TableCollection")

    def __iter__(self):
        '''Iterate over all the OTUs from all the tables

        Affected by the target_taxonomy instance variable, which narrows
        the scope of iteration to just those instances from that taxonomy
        self.target_taxonomy: list of str
            if not None, ignore those OTUs that are not from this clade,
            or more specific
        '''
        for table_types in (self.otu_table_objects, self.archive_table_objects):
            for table in table_types:
                for otu in table:
                    if self.target_taxonomy is None or otu.within_taxonomy(self.target_taxonomy):
                        yield otu

    def __len__(self):
        return sum(1 for e in self) # Why is this not automatic Python?

    def excluded_duplicate_distinct_genes(self):
        '''Filter the OTU table collection so that only a single OTU from each gene
        and sample combination is preserved, and iterate over the remaining
        OTUs.  When there is more than one different sequence for a different
        gene in the sample, then all of that gene from that sample are removed.
        Requires slurping the OTU table up.

        Returns
        -------
        A new OtuTableCollection object that has been filtered.

        '''
        sample_to_gene_to_otu = {}
        for otu in self:
            if otu.sample_name in sample_to_gene_to_otu:
                if otu.marker in sample_to_gene_to_otu[otu.sample_name]:
                    sample_to_gene_to_otu[otu.sample_name][otu.marker].append(otu)
                else:
                    sample_to_gene_to_otu[otu.sample_name][otu.marker] = [otu]
            else:
                sample_to_gene_to_otu[otu.sample_name] = {}
                sample_to_gene_to_otu[otu.sample_name][otu.marker] = [otu]

        for sample, gene_to_otu in sample_to_gene_to_otu.items():
            for gene, otus in gene_to_otu.items():
                logging.debug("Found %i OTUs for %s/%s" %(
                    len(otus), gene, otus[0].marker))
                if len(otus) == 1:
                    yield otus[0]

class StreamingOtuTableCollection:
    def __init__(self):
        self._otu_table_io_objects = []
        self._archive_table_io_objects = []
        self._otu_table_file_paths = []
        self._archive_table_file_paths = []

    def add_otu_table(self, input_otu_table_io):
        '''Add a regular style OTU table to the collection.

        Parameters
        ----------
        input_otu_table_ios: list of IO
            entries are open streams of OTU table data

        Returns
        -------
        None
        '''
        self._otu_table_io_objects.append(input_otu_table_io)

    def add_archive_otu_table(self, input_archive_table_io):
        self._archive_table_io_objects.append(input_archive_table_io)

    def add_otu_table_file(self, file_path):
        self._otu_table_file_paths.append(file_path)

    def add_archive_otu_table_file(self, file_path):
        self._archive_table_file_paths.append(file_path)

    def __iter__(self):
        '''Iterate over all the OTUs from all the tables. This can only be done once
        since the data is streamed in.
        '''
        for io in self._archive_table_io_objects:
            for otu in ArchiveOtuTable.read(io):
                yield otu
        for io in self._otu_table_io_objects:
            for otu in OtuTable.read(io):
                yield otu
        for file_path in self._archive_table_file_paths:
            for otu in ArchiveOtuTable.read(open(file_path)):
                yield otu
        for file_path in self._otu_table_file_paths:
            for otu in OtuTable.read(open(file_path)):
                yield otu
