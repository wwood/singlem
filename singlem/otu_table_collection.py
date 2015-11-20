from archive_otu_table import ArchiveOtuTable
from otu_table import OtuTable

class OtuTableCollection:
    def __init__(self):
        self.otu_table_objects = []
        self.archive_table_objects = []
        
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
        
    def __iter__(self, target_taxonomy=None):
        '''Iterate over all the OTUs from all the tables.
        
        Parameters
        ----------
        taxonomy: list of str
            if not None, ignore those OTUs that are not from this clade,
            or more specific
        '''
        for table_types in (self.otu_table_objects, self.archive_table_objects):
            for table in table_types:
                for otu in table:
                    if target_taxonomy is None or otu.within_taxonomy(target_taxonomy):
                        yield otu
