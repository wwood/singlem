import csv
import string
from archive_otu_table import ArchiveOtuTable

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
    def __init__(self):
        self.fields = string.split('gene sample sequence num_hits coverage taxonomy')
        self.data = []
        
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
            
    def write_to(self, output_io, fields_to_print):
        '''Output as a CSV file to the (open) I/O object
        
        Parameters
        ----------
        output_io: open io object
            this method neither opens nor closes this
        fields_to_print: list of str
            a list of names of fields to be printed
        '''
        if fields_to_print:
            field_indices_to_print = [self.fields.index(f) for f in fields_to_print]
        output_io.write("\t".join([self.fields[i] for i in field_indices_to_print])+"\n")            
        
        # When an element is actually multiple elements, join with a space,
        # and in any case convert everything to a string
        def to_printable(e):
            if hasattr(e, '__iter__'):
                return ' '.join([str(sub_e) for sub_e in e])
            elif isinstance(e, float):
                return "%.2f" % e
            else:
                return str(e)
        
        for d in self.data:
            output_io.write("\t".join([to_printable(d[i]) for i in field_indices_to_print])+"\n")
            
    def archive(self, singlem_packages):
        '''Return an archive object with the same data as this OTU table
        
        Parameters
        ----------
        singlem_packages: array of SingleMPackage objects
        
        Returns
        -------
        ArchiveOtuTable object
        '''
        archive = ArchiveOtuTable(singlem_packages)
        archive.fields = self.fields
        archive.data = self.data
        return archive

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