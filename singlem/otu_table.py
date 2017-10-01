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
        '''yield an OtuTableEntry object for each entry in the OTU table.
        This method is able to deal with streaming OTU tables.
        '''
        for i, d in enumerate(csv.reader(otu_table_io, delimiter="\t")):
            if len(d) < 5:
                raise Exception("Parse issue parsing line of OTU table: '%s'" % d)
            if i==0:
                fields = d
            else:
                if len(d) != len(fields):
                    raise Exception("Malformed OTU table detected, number of fields unexpected, on this line: %s" % str(d))
                e = OtuTableEntry()
                d[3] = int(d[3])
                d[4] = float(d[4])
                e.marker = d[0]
                e.sample_name = d[1]
                e.sequence = d[2]
                e.count = d[3]
                e.coverage = d[4]
                e.taxonomy = d[5]
                e.data = d
                e.fields = fields
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
            e.data = d
            e.fields = self.fields
            yield e

    def add(self, otu_table_entries):
        '''Add OtuTableEntry objects to this OTU table. Only the default
        data is saved'''
        for e in otu_table_entries:
            self.data.append([
                e.marker,
                e.sample_name,
                e.sequence,
                e.count,
                e.coverage,
                e.taxonomy])

    @staticmethod
    def read(input_otu_table_io):
        otus = OtuTable()
        for otu in OtuTable.each(input_otu_table_io):
            otus.data.append(otu.data)
        return otus

    def write_to(self, output_io, fields_to_print=DEFAULT_OUTPUT_FIELDS):
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
    def write_otus_to(otu_table_entries, output_io, fields_to_print=None):
        '''Output as a CSV file to the (open) I/O object

        Parameters
        ----------
        output_io: open io object
            this method neither opens nor closes this
        fields_to_print: list of str
            names of the fields to be printed. None indicates DEFAULT_OUTPUT_FIELDS
        '''

        if fields_to_print is None: fields_to_print = OtuTable.DEFAULT_OUTPUT_FIELDS
        output_io.write("\t".join(fields_to_print)+"\n")
        for o in otu_table_entries:
            field_indices_to_print = [o.fields.index(f) for f in fields_to_print]
            output_io.write("\t".join([OtuTable._to_printable(cell) for cell in [
                o.data[i] for i in field_indices_to_print]])+"\n")

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
