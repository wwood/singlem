import logging

class PlacementParser:
    def __init__(self, json, taxonomy_bihash, probability_threshold):
        self._json = json
        self._taxonomy_bihash = taxonomy_bihash
        self._sequence_to_placement = {}
        self._probability_threshold = probability_threshold
        for placement in json['placements']:
            for nm in placement['nm']:
                if nm[1] != 1:
                    raise Exception("Cannot handle jplace files with nm counts != 1")
                sequence = nm[0]
                if sequence in self._sequence_to_placement:
                    raise Exception(
                        "There appears to be duplicate names amongst placed sequences")
                self._sequence_to_placement[sequence] = placement

    def otu_placement(self, sequence_names):
        '''Return the most fully resolved taxonomy of the set of reads, pooling the
        placement confidences of each sequence. The returned taxonomy is the
        taxonomy whose placement probability sum is > len(sequence_names) *
        probability_threshold.

        Sequences not in the jplace are ignored. If a group is made up
        exclusively of sequences not in the jplace, None is returned.

        '''

        observed_parent_to_children = {}
        tax_probabilities = {}

        j = self._json
        classification_index = j['fields'].index('classification')
        likelihood_index = j['fields'].index('like_weight_ratio')
        root_tax = None

        # For each placement for each sequence
        for name in sequence_names:
            if name in self._sequence_to_placement:
                for p in self._sequence_to_placement[name]['p']:
                    # Get list of taxonomies from placed to root of tree
                    prob = p[likelihood_index]
                    tax = p[classification_index]
                    full_tax = []
                    while tax is not None:
                        full_tax.append(tax)
                        tax = self._taxonomy_bihash.child_to_parent[tax]
                    # Add that probability in the total hash
                    last_parent = None
                    for i, tax in enumerate(reversed(full_tax)):
                        if i == 0:
                            if root_tax is None:
                                root_tax = tax
                            elif tax != root_tax:
                                raise Exception(
                                    "Programming error - seem to have encountered 2 different roots")
                        if tax in tax_probabilities:
                            tax_probabilities[tax] += prob
                        else:
                            tax_probabilities[tax] = prob
                        if i != 0:
                            if last_parent in observed_parent_to_children:
                                if tax not in observed_parent_to_children[last_parent]:
                                    observed_parent_to_children[last_parent].append(tax)
                            else:
                                observed_parent_to_children[last_parent] = [tax]
                        last_parent = tax
            else:
                logging.debug("Skipping ORF %s as it does not seem to have been placed" % name)

        if root_tax is None:
            return None

        # Descend down the tree, picking the highest probability amongst the
        # children, or break if the probability is not greater than the
        # threshold.
        next_children = [root_tax]
        final_tax = []
        threshold = self._probability_threshold
        while True:
            # Find the max child and it's probability
            def get_prob(tax): return tax_probabilities[tax]
            max_child = max(next_children, key=get_prob)
            # Add it to the final tax if it is above the threshold
            if tax_probabilities[max_child] > threshold * len(sequence_names):
                final_tax.append(max_child)
                if max_child in observed_parent_to_children:
                    next_children = observed_parent_to_children[max_child]
                else:
                    break # If there is no children, break
            else:
                break # If below the threshold, break

        # Return the taxonomic placement
        return final_tax
