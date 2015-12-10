import csv

class UCFile:
    '''A parser for the UC format, or a subset thereof, defined at
    http://www.drive5.com/usearch/manual/ucout.html but meant particularly
    for parsing vsearch generated ones.
    '''
    
    def __init__(self, io):
        self.io = io
        
    def __iter__(self):
        '''
        Iterate over the UC input provided as an instance variable, 
        yielding UCEntry objects.
        '''
        for row in csv.reader(self.io, delimiter="\t"):
            if len(row) != 10:
                raise Exception("Unexpected format of UC file in this line: %s" % str(row))
            
            if row[0] != 'S': # don't yield the pontless S ones
                target = row[9]
                yield UCEntry(row[0], row[8], None if target=='*' else target)
                
class UCEntry:
    def __init__(self, record_type, query, target):
        self.record_type = record_type
        self.query = query
        self.target = target
