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
        return [self.marker, self.sample_name, self.sequence, self.count, self.coverage, self.taxonomy]

    def get_space_separated_field(self, field_name):
        '''
        Return a list of this field, which are separated by a space, and return them as a list
        '''
        try:
            field_index = self.fields.index(field_name)
        except ValueError:
            raise ValueError("Field named %s is not recorded in this OTU table entry" % field_name)
        return data[field_index].split(' ')
