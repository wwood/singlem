# Condense and taxonomic profiles

Marker genes produce separate views of a community. Condense combines their OTU evidence into one taxonomic profile so that marker-specific variation does not become duplicated community abundance.

The method considers markers together, applies trimmed summaries, and resolves abundance across the taxonomic hierarchy. This is why an OTU table is not itself a taxonomic profile and why profiles derived outside the main pipeline still need condensing.

Profile coverage is unfilled: a taxon's value excludes values assigned more specifically to its descendants. Filled coverage adds descendant values and is the basis for relative abundance.
