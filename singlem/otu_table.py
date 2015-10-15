import csv

class OtuTableEntry:
    marker = None
    sample_name = None
    sequence = None
    count = None
    taxonomy = None
    coverage = None
    
    def taxonomy_array(self):
        return self.taxonomy.split('; ')

class OtuTable:
    @staticmethod
    def each(otu_table_io):
        '''yield an OtuTableEntry object for each entry in the OTU table'''
        maybe_header = True
        for row in csv.reader(otu_table_io, delimiter="\t"):
            if row[0]=='gene' and maybe_header:
                maybe_header = False
                continue
            maybe_header = False
            
            if len(row) < 5:
                raise Exception("Parse issue parsing line of OTU table: '%s'" % row)
            e = OtuTableEntry()
            e.marker = row[0]
            e.sample_name = row[1]
            e.sequence = row[2]
            e.count = int(row[3])
            e.coverage = float(row[4])
            e.taxonomy = row[5]
            yield e

class TaxonomyTargetedOtuTable(OtuTable):
    def __init__(self, taxonomy):
        '''
        Parameters
        ----------
        taxonomy: list of str
            taxonomy, one entry in the list for each level
        '''
        self.target_taxonomy = taxonomy
        
    def each(self, otu_table_io):
        '''Like each(), except only yield those entries that belong to the
        given lineage.
        
        Parameters
        ----------
        otu_table_io: IO
            IO object of the OTU table
        '''
        for e in OtuTable().each(otu_table_io):
            if e.taxonomy_array()[:len(self.target_taxonomy)] == self.target_taxonomy:
                yield e