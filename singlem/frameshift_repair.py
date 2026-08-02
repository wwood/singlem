'''Repair frameshifts (single-base indels) in reads, using the frameshift
positions that DIAMOND already reports in its BTOP string during the prefilter.

Nanopore reads carry indels far more often than substitutions, and an indel is
much more costly to SingleM than a substitution: it shifts the reading frame, so
the translated (protein) alignment to the HMM breaks down and no window is
extracted at all. A substitution only changes one residue.

DIAMOND, when run with --frameshift, already finds these indels and records them
in the BTOP string as '/' and '\\' operations. This module turns those into an
edit of the nucleotide read that restores the reading frame:

  * '/' means the query is *missing* a base there, so a base is inserted. Which
    base was deleted is not recoverable from the read, so AMBIGUOUS_CHAR ('N') is
    inserted and resolved later by resolve_ambiguous_windows().
  * '\\' means the query has an *extra* base there, so that base is deleted.

BTOP coordinates
----------------
Positions are always in the query sequence's own coordinate frame, whichever
strand the alignment is on. On the reverse strand DIAMOND reports qstart > qend
and the walk proceeds downwards; the frameshift positions it yields still index
the query sequence as given, so the repair applies to it directly.

The per-operation query consumption is: a run of N matches consumes 3N
nucleotides; a mismatch consumes 3; a gap in the query ('-' as the first
character of the pair) consumes none; a gap in the subject consumes 3; and a
frameshift shifts the frame by one nucleotide, -1 for '/' and +1 for '\\'.
walk_btop() checks that this arithmetic lands exactly on qend, which is what
makes it safe to trust the positions.
'''

import logging
import re
from collections import Counter

# Base inserted where DIAMOND says the read is missing one. The identity of the
# deleted base is not recoverable from the read itself.
AMBIGUOUS_CHAR = 'N'

# Maximum number of mismatching (non-ambiguous) positions between an ambiguous
# window and a donor window for the donor to be usable.
DEFAULT_MAX_DIVERGENCE = 2

_BTOP_TOKEN = re.compile(r'(\d+)|([/\\])(.)|(.)(.)')


def parse_btop(btop):
    '''Yield the operations of a BTOP string in order.

    Yields tuples of:
      ('match', count)          count consecutive identical residues
      ('mismatch', query, subject)  a substitution or a gap on either side
      ('frameshift', char)      char is '/' (base missing from query) or '\\'
                                (extra base in query)
    '''
    position = 0
    while position < len(btop):
        match = _BTOP_TOKEN.match(btop, position)
        if match is None:
            raise ValueError(
                "Could not parse BTOP string '{}' at offset {}".format(btop, position))
        position = match.end()
        if match.group(1) is not None:
            yield ('match', int(match.group(1)))
        elif match.group(2) is not None:
            yield ('frameshift', match.group(2))
        else:
            yield ('mismatch', match.group(4), match.group(5))


def walk_btop(qstart, qend, btop):
    '''Walk a BTOP string and return the frameshifts it contains.

    Parameters
    ----------
    qstart, qend: int
        1-based query start/end as reported by DIAMOND. qstart > qend for
        alignments on the reverse strand.
    btop: str
        the BTOP string.

    Returns
    -------
    list of (position, kind), where position is a 0-based index into the query
    sequence and kind is 'deletion' (a base is missing from the read) or
    'insertion' (the read has an extra base). Returns None if the walk does not
    end exactly on qend, which means the BTOP semantics assumed here do not hold
    (e.g. a future DIAMOND change) and the read should be left alone.
    '''
    step = -1 if qstart > qend else 1
    position = qstart
    frameshifts = []
    try:
        operations = list(parse_btop(btop))
    except ValueError as e:
        logging.warning("%s", e)
        return None

    for operation in operations:
        if operation[0] == 'match':
            position += step * 3 * operation[1]
        elif operation[0] == 'frameshift':
            if operation[1] == '/':
                # On the reverse strand (step == -1), the missing base sits
                # immediately after `position` in the query's own (forward)
                # coordinates, not before it as on the forward strand.
                frameshifts.append((position if step == -1 else position - 1, 'deletion'))
                position -= step
            else:
                frameshifts.append((position - 1, 'insertion'))
                position += step
        else:
            # A gap in the query consumes no query nucleotides; everything else
            # (substitution, or a gap in the subject) consumes a codon.
            if operation[1] != '-':
                position += step * 3

    if position - step != qend:
        logging.warning(
            "BTOP string '%s' walked to query position %i but DIAMOND reported "
            "qend %i; not repairing frameshifts in this alignment",
            btop, position - step, qend)
        return None

    return frameshifts


def repair_frameshifts(sequence, frameshifts, ambiguous_char=AMBIGUOUS_CHAR):
    '''Apply frameshift repairs to a nucleotide sequence, restoring its reading
    frame.

    Deletions have ambiguous_char inserted (the deleted base is unknown);
    insertions have the extra base removed. Edits are applied from the 3' end
    backwards so that earlier positions stay valid.

    Parameters
    ----------
    sequence: str
        the nucleotide sequence to repair.
    frameshifts: list of (position, kind)
        as returned by walk_btop().

    Returns
    -------
    the repaired sequence.
    '''
    for (position, kind) in sorted(frameshifts, key=lambda f: -f[0]):
        if position < 0 or position > len(sequence):
            logging.warning(
                "Ignoring frameshift at out-of-range position %i of a %i bp sequence",
                position, len(sequence))
            continue
        if kind == 'deletion':
            sequence = sequence[:position] + ambiguous_char + sequence[position:]
        else:
            sequence = sequence[:position] + sequence[position + 1:]
    return sequence


def shift_position(position, frameshifts):
    '''Map a 0-based position in the original sequence to its position in the
    sequence returned by repair_frameshifts().'''
    shift = 0
    for (frameshift_position, kind) in frameshifts:
        if frameshift_position < position:
            shift += 1 if kind == 'deletion' else -1
    return position + shift


def _divergence(window, donor, ambiguous_char, max_mismatches=None):
    '''Number of mismatches between window and donor, ignoring positions where
    window is ambiguous.

    Returns None if the two are not comparable, or if max_mismatches is given and
    exceeded - counting stops as soon as that happens, since callers only care
    whether a donor is close enough to use.
    '''
    if len(window) != len(donor):
        return None
    mismatches = 0
    for (window_base, donor_base) in zip(window, donor):
        if window_base == ambiguous_char:
            continue
        if window_base != donor_base:
            mismatches += 1
            if max_mismatches is not None and mismatches > max_mismatches:
                return None
    return mismatches


# Blocks shorter than this match too many donors to narrow the search usefully,
# so no more blocks are cut than this allows however many ambiguous bases a
# window has.
_MIN_BLOCK_LENGTH = 4


class _DonorIndex:
    '''Donors indexed by exact-matching substring, so that an ambiguous window is
    compared only against plausible donors rather than all of them.

    If a window and a donor differ at no more than max_divergence unambiguous
    positions, and both are cut into max_divergence+1 blocks, at least one block
    must be free of mismatches (pigeonhole principle). A donor sharing no block
    with the window therefore cannot be within max_divergence, and can be skipped
    without being compared.

    Ambiguous bases complicate this: a block containing one can never match a
    donor exactly, so it is useless for lookup. Extra blocks are cut to
    compensate, enough that max_divergence+1 of them are expected to be
    ambiguity-free. A window with more ambiguous blocks than that allows for
    falls outside the guarantee, and candidates() then reports every donor rather
    than risk missing a valid one.
    '''

    def __init__(self, abundances, max_divergence, ambiguous_char, max_ambiguous):
        self._abundances = abundances
        self._max_divergence = max_divergence
        self._ambiguous_char = ambiguous_char

        donors_by_length = {}
        for donor in abundances:
            donors_by_length.setdefault(len(donor), []).append(donor)

        self._blocks_by_length = {}
        for (length, donors) in donors_by_length.items():
            num_blocks = min(max_divergence + 1 + max_ambiguous,
                             max(1, length // _MIN_BLOCK_LENGTH))
            bounds = [(i * length // num_blocks, (i + 1) * length // num_blocks)
                      for i in range(num_blocks)]
            blocks = [{} for _ in bounds]
            for donor in donors:
                for (block, (start, end)) in zip(blocks, bounds):
                    block.setdefault(donor[start:end], []).append(donor)
            self._blocks_by_length[length] = (bounds, blocks)

    def _by_descending_abundance(self, donors):
        # The sequence itself breaks ties, so the result does not depend on dict
        # ordering, matching the tie-break in resolve_ambiguous_windows().
        return sorted(donors, key=lambda donor: (-self._abundances[donor], donor))

    def candidates(self, window):
        '''Donors that could be within max_divergence of window, most abundant
        first. Donors of a different length are never comparable, so are
        omitted.'''
        entry = self._blocks_by_length.get(len(window))
        if entry is None:
            return ()
        (bounds, blocks) = entry

        candidates = {}
        num_usable_blocks = 0
        for (block, (start, end)) in zip(blocks, bounds):
            window_block = window[start:end]
            if self._ambiguous_char in window_block:
                continue
            num_usable_blocks += 1
            for donor in block.get(window_block, ()):
                candidates[donor] = None

        if num_usable_blocks <= self._max_divergence:
            # Too few unambiguous blocks for the pigeonhole guarantee to hold, so
            # a close donor might share none of them. Fall back to all donors.
            return self._by_descending_abundance(self._abundances)
        return self._by_descending_abundance(candidates)


def resolve_ambiguous_windows(window_sequences,
                              max_divergence=DEFAULT_MAX_DIVERGENCE,
                              ambiguous_char=AMBIGUOUS_CHAR,
                              repaired_names=None):
    '''Fill in ambiguous bases in window sequences from other windows in the same
    set.

    Repairing a deletion restores the reading frame but not the identity of the
    deleted base, so the window carries an ambiguous_char. Such a window cannot
    match a reference window exactly, and would be counted as its own OTU. Here
    each ambiguous window takes its missing base(s) from the most abundant
    unambiguous window that is within max_divergence mismatches of it, which is
    almost always the same organism's window recovered from reads that happened
    not to have an indel there.

    Windows with no such neighbour are left as they are: an ambiguous window
    still contributes coverage, and inventing a base for it would be worse than
    admitting the uncertainty.

    Parameters
    ----------
    window_sequences: list of UnalignedAlignedNucleotideSequence
        modified in place (their aligned_sequence attribute is updated).
    max_divergence: int
        maximum number of mismatches at unambiguous positions.
    repaired_names: set of str, or None
        if given, only windows whose .name is in this set are eligible to be
        resolved, since an ambiguous_char elsewhere is not from frameshift
        repair (e.g. an 'N' already present in the raw read) and inventing a
        base for it would misrepresent the read. If None, every ambiguous
        window is eligible, as when no such provenance is available.

    Returns
    -------
    the number of window sequences that were resolved.
    '''
    ambiguous = [s for s in window_sequences
                 if ambiguous_char in s.aligned_sequence
                 and (repaired_names is None or s.name in repaired_names)]
    if len(ambiguous) == 0:
        return 0

    abundances = Counter(
        s.aligned_sequence for s in window_sequences
        if ambiguous_char not in s.aligned_sequence)
    if len(abundances) == 0:
        logging.debug(
            "No unambiguous window sequences available to resolve %i ambiguous one(s)",
            len(ambiguous))
        return 0

    index = _DonorIndex(
        abundances, max_divergence, ambiguous_char,
        max_ambiguous=max(s.aligned_sequence.count(ambiguous_char) for s in ambiguous))

    # Cache by window sequence, since many reads usually share one.
    resolutions = {}
    num_resolved = 0
    for sequence in ambiguous:
        window = sequence.aligned_sequence
        if window not in resolutions:
            best_donor = None
            best_key = None
            for donor in index.candidates(window):
                abundance = abundances[donor]
                # Candidates come most abundant first and abundance dominates the
                # key below, so once a donor is less abundant than the best found
                # no later one can beat it.
                if best_key is not None and -abundance > best_key[0]:
                    break
                divergence = _divergence(
                    window, donor, ambiguous_char, max_mismatches=max_divergence)
                if divergence is None:
                    continue
                # Most abundant wins; the sequence itself breaks ties so the
                # result does not depend on dict ordering.
                key = (-abundance, divergence, donor)
                if best_key is None or key < best_key:
                    best_key = key
                    best_donor = donor
            if best_donor is None:
                resolutions[window] = None
            else:
                resolutions[window] = ''.join(
                    donor_base if window_base == ambiguous_char else window_base
                    for (window_base, donor_base) in zip(window, best_donor))
                logging.debug(
                    "Resolved ambiguous window %s to %s using a donor of abundance %i",
                    window, resolutions[window], -best_key[0])
        resolved = resolutions[window]
        if resolved is not None:
            sequence.aligned_sequence = resolved
            num_resolved += 1

    return num_resolved
