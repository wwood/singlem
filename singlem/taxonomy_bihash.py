import csv

class TaxonomyBihash:
    def __init__(self):
        self.parent_to_children = {}
        self.child_to_parent = {}

    @staticmethod
    def parse_taxtastic_taxonomy(taxtastic_taxonomy_io):
        '''Read a taxtastic style taxonomy in, given an open IO object to e.g. the
        file'''

        first = True
        csv_reader = csv.reader(taxtastic_taxonomy_io, delimiter=',')
        bihash = TaxonomyBihash()
        parent_to_children = bihash.parent_to_children
        child_to_parent = bihash.child_to_parent
        for row in csv_reader:
            if first:
                first = False
                continue
            tax_id = row[0]
            parent_id = row[1]

            if tax_id == 'Root':
                child_to_parent[tax_id] = None
            else:
                if parent_id in parent_to_children:
                    parent_to_children[parent_id].append(tax_id)
                else:
                    parent_to_children[parent_id] = [tax_id]
                if tax_id in child_to_parent:
                    raise Exception(
                        "Found duplicate parents for child ID %s when parsing taxonomy file" %
                        tax_id)
                child_to_parent[tax_id] = parent_id
        return bihash
