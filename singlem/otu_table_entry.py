from .taxonomy import TaxonomyUtils

class OtuTableEntry:
    marker = None
    sample_name = None
    sequence = None
    count = None
    taxonomy = None
    coverage = None
    data = None
    fields = None

    def taxonomy_array(self):
        return TaxonomyUtils.split_taxonomy(self.taxonomy)

    def within_taxonomy(self, target_taxonomy):
        '''Return true iff the OTU has been assigned within this taxonomy,
        else false

        Parameters
        ----------
        taxonomy: list of str
            each taxonomy level
        '''
        return (self.taxonomy_array()[:len(target_taxonomy)] == target_taxonomy)

    def add_found_data(self, found_in):
        if 'found_in' not in self.fields:
            self.fields = self.fields.copy()  # Create a copy to avoid global modifications
            self.fields.append('found_in')
            self.data.append(found_in)
        else:
            try:
                self.data[self.fields.index('found_in')] = self.data[self.fields.index('found_in')] + ',' + found_in
            except IndexError:
                self.data.append(found_in)

    def __str__(self):
        return "\t".join([self.marker, self.sample_name, self.sequence,
                          str(self.count), str(self.coverage), self.taxonomy])

    def to_list(self):
        return [self.marker, self.sample_name, self.sequence,
            self.count, self.coverage, self.taxonomy]
