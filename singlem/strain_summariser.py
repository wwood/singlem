from otu_table import OtuTableEntry


class DifferenceOTUEntry(OtuTableEntry):
    difference_in_bp = None
    
    @staticmethod
    def create_from_otu_table_entry(otu_table_entry):
        '''Create a new entry out of an OtuTableEntry'''
        d = DifferenceOTUEntry()
        d.marker = otu_table_entry.marker
        d.sample_name = otu_table_entry.sample_name
        d.sequence = otu_table_entry.sequence
        d.count = otu_table_entry.count
        d.taxonomy = otu_table_entry.taxonomy
        d.coverage = otu_table_entry.coverage
        return d

class StrainSummariser:
    def summarise_strains(self, **kwargs):
        table_collection = kwargs.pop('table_collection')
        output_table_io = kwargs.pop('output_table_io')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        output_table_io.write("\t".join(['type',
                                         'gene',
                                         'sample',
                                         'difference_in_bp',
                                         'sequence',
                                         'num_hits',
                                         'coverage',
                                         'taxonomy'
                                         ])+"\n")
        last_sample_name = None
        last_gene = None
        current_sample_entries = []
        
        for e in table_collection:
            if last_sample_name is None:
                last_sample_name = e.sample_name
                last_gene = e.marker
                current_sample_entries.append(e)
            elif e.sample_name == last_sample_name and e.marker == last_gene:
                current_sample_entries.append(e)
            else:
                self._process_sample(current_sample_entries, output_table_io)
                last_sample_name = e.sample_name
                last_gene = e.marker
                current_sample_entries = [e]
        self._process_sample(current_sample_entries, output_table_io)
                 
    def _process_sample(self, entries, output_table_io):
        # pick the reference sequence based on abundance - most abundant
        # first, secondly by alphabetical order of the sequence
        reference_entry = max(entries, key=lambda e: [e.coverage, e.taxonomy])
        output_table_io.write("\t".join(['reference',
                                         reference_entry.marker,
                                         reference_entry.sample_name,
                                         '0',
                                         reference_entry.sequence,
                                         str(reference_entry.count),
                                         str(reference_entry.coverage),
                                         reference_entry.taxonomy
                                         ])+"\n")
        for e in sorted(self._differences(reference_entry, entries),
                        key=lambda d: -d.difference_in_bp):
            output_table_io.write("\t".join([
                                         'strain',
                                         e.marker,
                                         e.sample_name,
                                         str(e.difference_in_bp),
                                         e.sequence,
                                         str(e.count),
                                         str(e.coverage),
                                         e.taxonomy    
                                         ])+"\n")
        
    def _differences(self, reference_otu, otu_table_iterator):
        '''Take an OTU table and an entry from it. Return a list of DifferenceOTUEntry
        entry that includes the difference between the OTU and the reference
        entry in percent identity.
        
        Ignores the reference OTU if it is present in the otu_table
        
        Parameters
        ----------
        otu_table: OtuTable
            OTUs to compare to the reference
        reference_otu: OtuTableEntry
            the reference OTU
            
        Returns
        -------
        generator over DifferenceOTUEntry objects
        '''
        for otu in otu_table_iterator:
            # ignore the reference ones
            if otu.sequence == reference_otu.sequence: continue
            
            d = DifferenceOTUEntry.create_from_otu_table_entry(otu)
            d.difference_in_bp = 0
            for i, ref in enumerate(reference_otu.sequence):
                if otu.sequence[i] != ref:
                    d.difference_in_bp += 1
            yield d
