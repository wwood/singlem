This mode is an integral part of creating a new SingleM package from scratch. The purpose of the mode is to choose where in the HMM is the best place for a starting window. At the moment this means the most conserved stretch.

The best position output is the one that has the most nucleic acids that overlap the HMM at the window starting at that position.

The input is an alignment created by applying the HMMER tool `hmmalign`, which is converted to FASTA format using [seqmagick](https://github.com/fhcrc/seqmagick) `convert`.

The window positions in the default SingleM packages were chosen through this method, supplying an alignment created from the reads of a complex soil metagenome. Whole gene sequences may also be appropriate as input (after alignment).

Once a best window position is chosen through this process, [singlem create](/advanced/create) is used to finalise creation of the SingleM package.
