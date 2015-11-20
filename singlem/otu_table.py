import csv
import string
from archive_otu_table import ArchiveOtuTable
from otu_table_entry import OtuTableEntry

class OtuTable:
    DEFAULT_OUTPUT_FIELDS = string.split('gene sample sequence num_hits coverage taxonomy')
    
    def __init__(self):
        self.fields = self.DEFAULT_OUTPUT_FIELDS
        self.data = []
        
    @staticmethod
    def each(otu_table_io):
        '''yield an OtuTableEntry object for each entry in the OTU table'''
        otus = OtuTable.read(otu_table_io)
        for e in otus:
            yield e
            
    def __iter__(self):
        for d in self.data:
            e = OtuTableEntry()
            e.marker = d[0]
            e.sample_name = d[1]
            e.sequence = d[2]
            e.count = d[3]
            e.coverage = d[4]
            e.taxonomy = d[5]
            yield e
            
    @staticmethod
    def read(input_otu_table_io):
        otus = OtuTable()
        for i, row in enumerate(csv.reader(input_otu_table_io, delimiter="\t")):
            if len(row) < 5:
                raise Exception("Parse issue parsing line of OTU table: '%s'" % row)
            if i==0:
                otus.fields = row
            else:
                if len(row) != len(otus.fields):
                    raise Exception("Malformed OTU table detected, number of fields unexpected, on this line: %s" % str(row))
                row[3] = int(row[3])
                row[4] = float(row[4])
                otus.data.append(row)
        return otus
            
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
        
        for d in self.data:
            output_io.write("\t".join([self._to_printable(d[i]) for i in field_indices_to_print])+"\n")
            
    @staticmethod
    def write_otus_to(otu_table_entries, output_io):
        '''Output as a CSV file to the (open) I/O object
        
        Parameters
        ----------
        output_io: open io object
            this method neither opens nor closes this
        '''
        output_io.write("\t".join(OtuTable.DEFAULT_OUTPUT_FIELDS)+"\n")
        
        for d in otu_table_entries:
            output_io.write("\t".join([OtuTable._to_printable(cell) for cell in [\
                d.marker,
                d.sample_name,
                d.sequence,
                d.count,
                d.coverage,
                d.taxonomy]])+"\n")
         
    @staticmethod
    def _to_printable(e):
        '''For printing OtuTableEntry parts
        When an element is actually multiple elements, join with a space,
        and in any case convert everything to a string'''
        if hasattr(e, '__iter__'):
            return ' '.join([str(sub_e) for sub_e in e])
        elif isinstance(e, float):
            return "%.2f" % e
        else:
            return str(e)
            
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
