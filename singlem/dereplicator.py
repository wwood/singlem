import random

from graftm.greengenes_taxonomy import GreenGenesTaxonomy

class Dereplicator:
    '''Return a cut down list of identifiers, where only a single representative
    is chosen from some level of the taxonomy.

    Params
    ------
    input_identifiers: list of str
        Identifiers to dereplicate
    dereplication_level: int
        Level for dereplication
    taxonomy_hash: str to list of str
        Id to taxonomy as list
    preferred_lineages: list of str
        Given the choice for a dereplication group, choose one of these
        lineages as the representative.

    '''
    def dereplicate(self, input_identifiers, dereplication_level,
                    taxonomy_hash, preferred_lineages):
        to_return = []
        dereplication_groups = {}
        for identifier in input_identifiers:
            tax = taxonomy_hash[identifier]
            if len(tax) < dereplication_level:
                to_return.append(identifier) # Do not dereplicate when there's
                                             # insufficient taxonomy.
            else:
                entry = '; '.join(tax[:(dereplication_level+1)])
                if entry in dereplication_groups:
                    dereplication_groups[entry].append(identifier)
                else:
                    dereplication_groups[entry] = [identifier]

        preferred_set = set(preferred_lineages)
        for group in dereplication_groups.values():
            preferreds = list([i for i in group if i in preferred_set])
            if len(preferreds) > 0:
                to_return.append(random.choice(preferreds))
            else:
                to_return.append(random.choice(group))

        return to_return
