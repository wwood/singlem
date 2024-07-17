import logging

from .archive_otu_table import ArchiveOtuTable
from .otu_table import OtuTable
from .utils import FastaNameToSampleName


class GenomeFastaMux:
    def __init__(self, genome_fasta_files):
        # Create dict from enumerate of strings
        self._genomes = []
        self._fasta_to_index = {}
        for i, genome_fasta_file in enumerate(genome_fasta_files):
            self._genomes.append(genome_fasta_file)
            if genome_fasta_file in self._fasta_to_index:
                raise ValueError("Duplicate genome fasta file: %s" % genome_fasta_file)
            self._fasta_to_index[genome_fasta_file] = i

    def fasta_to_prefix(self, fasta_file):
        return self._fasta_to_index[fasta_file]

    def demux_otu_table(self, otu_table_object):
        logging.info("Demultiplexing genomes from OTU table")
        archive = otu_table_object.archive(None)
        new_otus = []

        if ArchiveOtuTable.version != 4:
            raise ValueError("Unsupported OTU table version provided to genome demux - likely a programming bug: %d" % ArchiveOtuTable.version)

        for otu in archive:
            # Add an OTU for each genome
            genome_to_otus = {}
            for (sequence_name,
                    nucleotides_aligned,
                    read_unaligned_sequence) in zip(
                    otu.read_names(),
                    otu.nucleotides_aligned(),
                    otu.read_unaligned_sequences()
                    ):
                genome_id = int(sequence_name.split('|')[0])
                original_sequence_name = '|'.join(sequence_name.split('|')[1:])
                genome_name = FastaNameToSampleName.fasta_to_name(self._genomes[genome_id])
                if genome_name not in genome_to_otus:
                    genome_to_otus[genome_name] = []
                genome_to_otus[genome_name].append([
                    otu.marker,
                    genome_name,
                    otu.sequence,
                    1,
                    1.0,  # Hard to calculate coverage, and it is meaningless, so eh.
                    otu.taxonomy,
                    [original_sequence_name],
                    [nucleotides_aligned],
                    otu.taxonomy_by_known(),
                    [read_unaligned_sequence],
                    otu.equal_best_hit_taxonomies(),
                    otu.taxonomy_assignment_method()
                ])
            for genome_name, otus in genome_to_otus.items():
                genome_otu = [
                    otus[0][0],
                    otus[0][1],
                    otus[0][2],
                    sum([otu[3] for otu in otus]),
                    sum([otu[4] for otu in otus]),
                    otus[0][5],
                    [read_name for otu in otus for read_name in otu[6]],
                    [nucleotides_aligned for otu in otus for nucleotides_aligned in otu[7]],
                    otus[0][8],
                    [read_unaligned_sequence for otu in otus for read_unaligned_sequence in otu[9]],
                    otus[0][10],
                    otus[0][11]
                ]
                new_otus.append(genome_otu)

        # Sort by sample then gene, because standard not to intermix samples
        new_otus = sorted(new_otus, key=lambda otu: (otu[1], otu[0]))
        to_return = OtuTable()
        to_return.data = new_otus
        return to_return
