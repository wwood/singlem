PPLACER_ASSIGNMENT_METHOD = 'pplacer'
DIAMOND_ASSIGNMENT_METHOD = 'diamond'
DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD = 'diamond_example'
ANNOY_ASSIGNMENT_METHOD = 'annoy'
ANNOY_THEN_DIAMOND_ASSIGNMENT_METHOD = 'annoy_then_diamond'
SCANN_THEN_DIAMOND_ASSIGNMENT_METHOD = 'scann_then_diamond'
SCANN_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD = 'scann_naive_then_diamond'
SMAFA_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD = 'smafa_naive_then_diamond'
NO_ASSIGNMENT_METHOD = 'no_assign_taxonomy'

QUERY_BASED_ASSIGNMENT_METHOD = 'singlem_query_based'

class TaxonomyUtils:
    @staticmethod
    def split_taxonomy(taxonomy_string):
        if taxonomy_string:
            tax = [t.strip() for t in taxonomy_string.split(';')]
            while len(tax) > 0 and tax[-1] == '':
                tax = tax[:-1]
            return tax
        else:
            return None

    @staticmethod
    def lca_taxonomy_of_strings(taxonomy_strings):
        hit_taxonomies = list([list([ta.strip() for ta in t.split(';') if ta.strip() != '']) for t in taxonomy_strings])
        return TaxonomyUtils.lca_taxonomy_of_taxon_lists(hit_taxonomies)

    @staticmethod
    def lca_taxonomy_of_taxon_lists(taxonomy_lists):
        lca = taxonomy_lists[0]
        for taxonomy in taxonomy_lists[1:]:
            if len(taxonomy) < len(lca):
                lca = lca[:len(taxonomy)]
            for i, tax in enumerate(taxonomy):
                if i >= len(lca) or tax != lca[i]:
                    lca = lca[:i]
                    break
        return '; '.join(lca)
