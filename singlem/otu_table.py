import csv

from .archive_otu_table import ArchiveOtuTable
from .otu_table_entry import OtuTableEntry


class OtuTable:
    DEFAULT_OUTPUT_FIELDS = str.split('gene sample sequence num_hits coverage taxonomy')

    def __init__(self):
        self.fields = self.DEFAULT_OUTPUT_FIELDS
        self.data = []

    @classmethod
    def _clear_cache(cls):
        # For testing only, to clear the class-variable cache
        cls.DEFAULT_OUTPUT_FIELDS = str.split('gene sample sequence num_hits coverage taxonomy')

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
                try:
                    d[3] = int(d[3])
                except ValueError:
                    raise Exception("Malformed OTU table detected, num_hits column is not an integer, on line %i: %s" % (i+1, str(d)))
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

    def add_with_extras(self, otu_table_entries, extra_entries):
        '''Add OtuTableEntry objects to this OTU table. Default data
        and extra entries are saved'''
        for e in otu_table_entries:
            e_list = [
                e.marker,
                e.sample_name,
                e.sequence,
                e.count,
                e.coverage,
                e.taxonomy]
            for extra_field in extra_entries:
                try:
                    e_list.append(e.data[e.fields.index(extra_field)])
                except IndexError:
                    e_list.append('')
            self.data.append(e_list)

        for extra_field in extra_entries:
            if extra_field not in self.fields:
                self.fields.append(extra_field)

    def add_extras_no_data(self, extra_entries):
        '''Add extra entries to empty OTU table.'''
        if self.data:
            raise Exception("Cannot add extra entries to an OTU table that already has data")

        for extra_field in extra_entries:
            if extra_field not in self.fields:
                self.fields.append(extra_field)

    def sort_by_marker(self):
        '''Sort the OTU table by marker gene.
        '''
        marker_column = self.fields.index('gene')
        self.data.sort(key=lambda x: x[marker_column])

    @staticmethod
    def read(input_otu_table_io):
        otus = OtuTable()
        for otu in OtuTable.each(input_otu_table_io):
            otus.data.append(otu.data)
        return otus

    def write_to(self, output_io, fields_to_print=DEFAULT_OUTPUT_FIELDS, print_header=True):
        '''Output as a CSV file to the (open) I/O object

        Parameters
        ----------
        output_io: open io object
            this method neither opens nor closes this
        fields_to_print: list of str
            a list of names of fields to be printed
        '''
        field_indices_to_print = [self.fields.index(f) for f in fields_to_print]
        if print_header:
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
        if hasattr(e, '__iter__') and not isinstance(e, str):
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

    def rename_samples(self, renaming_dict):
        '''Rename samples according to the hash, mutating the state of the data in this object
        '''
        sample_column = self.fields.index('sample')
        for d in self.data:
            if d[sample_column] in renaming_dict:
                d[sample_column] = renaming_dict[d[sample_column]]