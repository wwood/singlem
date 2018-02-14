from taxonomy import Taxonomy

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
        return Taxonomy.split_taxonomy(self.taxonomy)

    def within_taxonomy(self, target_taxonomy):
        '''Return true iff the OTU has been assigned within this taxonomy,
        else false

        Parameters
        ----------
        taxonomy: list of str
            each taxonomy level
        '''
        return (self.taxonomy_array()[:len(target_taxonomy)] == target_taxonomy)

    def __str__(self):
        return "\t".join([self.marker, self.sample_name, self.sequence,
                          str(self.count), str(self.coverage), self.taxonomy])
