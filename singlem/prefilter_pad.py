import logging

from .sequence_classes import SeqReader

# Number of residues of flanking sequence to include on each side of the window.
DEFAULT_FLANK_LENGTH = 30
# Character used to pad sequences that do not extend far enough, and to fill
# deletions (gaps) within the window.
PAD_CHAR = 'X'


def window_alignment_positions(aligned_sequences, window_start_position, window_num_positions):
    '''Return the list of alignment column indices that make up the window, given
    a set of sequences aligned to the same HMM.

    Insert columns are determined globally across all sequences (a column is an
    insert column if any sequence has a lower-case residue there), matching the
    behaviour of MetagenomeOtuFinder so that the window here is identical to the
    one used by 'pipe'. This is necessary because hmmalign output represents both
    match-state deletions and insert-state gaps as '-', so match columns cannot be
    identified from a single sequence.

    Parameters
    ----------
    aligned_sequences: list of Sequence / AlignedProteinSequence
        sequences aligned to the HMM (all the same length).
    window_start_position: int
        0-based index of the first window match state, in match-state coordinates
        (i.e. not counting insert columns). This is the package's
        singlem_position().
    window_num_positions: int
        number of match states (amino acids) in the window, i.e. window_size() / 3.

    Returns
    -------
    list of int of length window_num_positions: the alignment column indices of
    the window.
    '''
    from .metagenome_otu_finder import MetagenomeOtuFinder

    finder = MetagenomeOtuFinder()
    ignored_columns = finder._find_lower_case_columns(aligned_sequences)
    start_column = finder._upper_case_position_to_alignment_position(
        window_start_position, ignored_columns)
    return finder._best_position_to_chosen_positions(
        start_column, window_num_positions, ignored_columns)


def pad_aligned_sequence(aligned_sequence, window_columns,
                         flank_length=DEFAULT_FLANK_LENGTH, pad_char=PAD_CHAR):
    '''Convert an HMM-aligned protein sequence into a fixed-length padded
    sequence of length (flank_length + len(window_columns) + flank_length), with
    the window occupying a consistent set of positions in the output.

    The output is composed of three parts:
      * flank_length residues immediately N-terminal to the window (inserts
        included), left-padded with pad_char if the sequence is too short;
      * the window itself: exactly len(window_columns) residues taken from the
        window's alignment columns, with inserts within the window removed
        (because they fall in columns not listed in window_columns) and any
        deletions replaced by pad_char;
      * flank_length residues immediately C-terminal to the window (inserts
        included), right-padded with pad_char if the sequence is too short.

    All residues are upper-cased in the output.

    Parameters
    ----------
    aligned_sequence: str
        the aligned sequence (a single row from hmmalign output). Residues are
        upper- or lower-case letters, gaps are '-' or '.'.
    window_columns: list of int
        ascending alignment column indices of the window, as returned by
        window_alignment_positions().
    flank_length: int
        number of flanking residues to include on each side.
    pad_char: str
        single character used for padding and to fill deletions in the window.

    Returns
    -------
    str of length (2*flank_length + len(window_columns)), or None if the sequence
    does not span the window (i.e. the first or last window column is a gap).
    '''
    first_window_column = window_columns[0]
    last_window_column = window_columns[-1]

    # Require the sequence to actually cover both ends of the window, matching the
    # behaviour of MetagenomeOtuFinder.find_windowed_sequences.
    if aligned_sequence[first_window_column] in '-.' or \
            aligned_sequence[last_window_column] in '-.':
        return None

    window = ''.join(
        pad_char if aligned_sequence[i] in '-.' else aligned_sequence[i].upper()
        for i in window_columns)

    left_residues = [c for c in aligned_sequence[:first_window_column] if c.isalpha()]
    right_residues = [c for c in aligned_sequence[last_window_column + 1:] if c.isalpha()]

    left_flank = ''.join(left_residues[-flank_length:]).upper().rjust(flank_length, pad_char)
    right_flank = ''.join(right_residues[:flank_length]).upper().ljust(flank_length, pad_char)

    return left_flank + window + right_flank


class PrefilterPadder:
    '''Produce a prefilter FASTA in which every on-target sequence is padded to a
    fixed length, with the SingleM window occupying a consistent set of
    positions. This allows DIAMOND to be configured to reject alignments that do
    not cover the window.'''

    def pad_metapackage(self, **kwargs):
        metapackage = kwargs.pop('metapackage')
        output_fasta = kwargs.pop('output_fasta')
        flank_length = kwargs.pop('flank_length', DEFAULT_FLANK_LENGTH)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments: {}".format(kwargs))

        # Only sequences that are actually in the prefilter DIAMOND database
        # should be padded, so the padded FASTA corresponds exactly to the
        # (CD-HIT dereplicated) prefilter DB rather than the full set of
        # on-target sequences.
        prefilter_names = self._prefilter_sequence_names(metapackage)

        total_written_seqs = 0
        with open(output_fasta, 'w') as out:
            for pkg in metapackage:
                if not pkg.is_protein_package():
                    logging.debug(
                        "Skipping non-protein package {}".format(pkg.base_directory()))
                    continue
                total_written_seqs += self._pad_package(
                    pkg, out, flank_length, prefilter_names)
        logging.info("Wrote {} padded sequences in total".format(total_written_seqs))

    def _prefilter_sequence_names(self, metapackage):
        '''Return the set of sequence names present in the metapackage's prefilter
        DIAMOND database.'''
        import extern

        prefilter_db_path = metapackage.prefilter_db_path()
        if prefilter_db_path is None:
            raise Exception(
                "The metapackage has no prefilter DIAMOND database, so "
                "prefilter-pad cannot restrict output to it")

        logging.info("Reading sequence names from prefilter DB {} ..".format(
            prefilter_db_path))
        names = set()
        stdout = extern.run("diamond getseq --db '{}'".format(prefilter_db_path))
        for line in stdout.splitlines():
            if line.startswith('>'):
                names.add(line[1:].split()[0])
        logging.info("Found {} sequences in the prefilter DB".format(len(names)))
        return names

    def _pad_package(self, pkg, out, flank_length, prefilter_names):
        # Import here to avoid the heavy pipe import chain when this module is
        # imported for the standalone padding functions only.
        from .pipe_sequence_extractor import _align_proteins_to_hmm

        if pkg.version < 3:
            raise Exception(
                "Padding a prefilter DB only works on version 3+ SingleM packages")

        window_start_position = pkg.singlem_position()
        window_num_positions = int(pkg.window_size() / 3)

        logging.info("Reading FASTA from {} ..".format(pkg.base_directory()))

        # Select only sequences that are present in the prefilter DB. Membership
        # in prefilter_names already implies the sequence is on-target and passed
        # the prefilter's own X-containing filter, so no further filtering by
        # target domain is needed here.
        sequences = []
        with open(pkg.graftm_package().unaligned_sequence_database_path()) as f:
            for (name, seq, _) in SeqReader().readfq(f):
                if name in prefilter_names:
                    sequences.append((name, seq))
        logging.info("Found {} sequences also present in the prefilter DB".format(
            len(sequences)))

        if len(sequences) == 0:
            return 0

        aligned = _align_proteins_to_hmm(
            sequences, pkg.graftm_package().alignment_hmm_path())
        if len(aligned) == 0:
            return 0

        window_columns = window_alignment_positions(
            aligned, window_start_position, window_num_positions)

        written = 0
        for aligned_sequence in aligned:
            padded = pad_aligned_sequence(
                aligned_sequence.seq, window_columns, flank_length=flank_length)
            if padded is None:
                logging.debug(
                    "Sequence {} does not span the window, skipping".format(
                        aligned_sequence.name))
                continue
            out.write(">{}\n{}\n".format(aligned_sequence.name, padded))
            written += 1
        logging.info("Wrote {} padded sequences for this package".format(written))
        return written
