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
        
    def set_target_taxonomy_by_string(self, taxonomy_string):
        '''Set the target_taxonomy instance variable by a string, which
        gets parsed into the requisite array form and stored in the instance
        variable'''
        self.target_taxonomy = Taxonomy.split_taxonomy(taxonomy_string)
        
    def __iter__(self):
        '''Iterate over all the OTUs from all the tables.
        
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
