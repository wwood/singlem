class Taxonomy:
    @staticmethod
    def split_taxonomy(taxonomy_string):
        if taxonomy_string:
            tax = [t.strip() for t in taxonomy_string.split(';')]
            while len(tax) > 0 and tax[-1] == '':
                tax = tax[:-1]
            return tax
        else:
            return None
