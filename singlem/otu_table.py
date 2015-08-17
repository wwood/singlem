import csv

class OtuTableEntry:
    marker = None
    sample_name = None
    sequence = None
    count = None
    taxonomy = None
    

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
            
            if len(row) < 5: raise Exception("Parse issue parsing line of OTU table: '%s'" % row)
            e = OtuTableEntry()
            e.marker = row[0]
            e.sample_name = row[1]
            e.sequence = row[2]
            e.count = int(row[3])
            e.taxonomy = row[4]
            yield e
