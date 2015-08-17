import otu_table
import logging

class KnownOtuTable:
    def parse_otu_tables(self, otu_table_list):
        '''Given a list of paths to OTU tables, parse them in to this object.
        If the same OTU sequence is encountered multiple times, the first one
        encountered is trusted. So the list provided should be in descending
        order of utility / trustworthiness
        
        otu_table_list: list of str
            paths to the OTU tables'''
        self._sequence_to_table_entry = {}
        for table_path in otu_table_list:
            count = 0
            for entry in otu_table.OtuTable.each(open(table_path)):
                key = entry.sequence
                if key in self._sequence_to_table_entry:
                    logging.debug("Ignoring %s when parsing in the OTU table %s as it has been seen previously" % (entry.sequence,
                                                                                                                   table_path))
                else:
                    self._sequence_to_table_entry[key] = entry
                    count += 1
            logging.debug("Parsed in %i entries from OTU table %s" % (count, table_path))
            
    def __getitem__(self, sequence):
        return self._sequence_to_table_entry[sequence]
    
    def __len__(self):
        return len(self._sequence_to_table_entry)
    
    def __contains__(self, sequence):
        return sequence in self._sequence_to_table_entry
