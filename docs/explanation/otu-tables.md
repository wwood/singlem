# OTU tables

An operational taxonomic unit table records the short marker-window sequences found in each sample. Each row identifies the marker, sample, nucleotide sequence, number of supporting reads, estimated genome coverage, and an initial taxonomic assignment.

OTU sequences exist independently of taxonomy. This makes them useful for comparing communities even when reference taxonomies change or a sequence belongs to a novel lineage.

An extended OTU table also retains read names, alignment details, and unaligned sequence. An archive OTU table stores the richer information in a versioned machine-readable form. Archive tables support later reannotation and database construction.

Coverage estimates genome-wide sequence coverage from reads observed across a marker window. It is distinct from the raw number of matching reads.
