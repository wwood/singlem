'''The archive OTU table: the richer JSON form of an OTU table, which retains
per-read data (names, unaligned sequences, aligned lengths) so that `renew` can
re-run from it.

On-disk encodings
-----------------
Versions 1-4 store `otus` as a list of rows, each row a list of values in
`fields` order. Version 5 stores the same rows column-wise, hoists fields whose
value is identical in every row into `constant_fields`, and dereplicates read
sequences into a shared `reads` list that the rows refer to. Those three changes
only remove repetition; nothing is added or lost.

Version 5 also adds one field, `reads_with_repaired_deletions`, which records
which of an OTU's reads had an ambiguous base inserted into their window by
frameshift repair. `pipe` needs that to tell such a base apart from one already
present in a raw read, and until version 5 it only existed in memory, so the
resolution of ambiguous windows had to happen inside the same `pipe` run. It is
None whenever no read of the OTU has one, which is every OTU of a run without
frameshift repair, so constant hoisting collapses the whole column to a single
entry there.

The gain is real but modest: on test/data/SRR8653040.json (95 OTUs, 1197 reads)
gzip -9 goes from 36.6 KB to 35.6 KB and zstd -19 from 32.1 KB to 30.4 KB. Most
of what remains is beyond the reach of any layout - the read sequences
themselves, and the per-package sha256s in the header, which are incompressible
and alone account for 13% of that gzipped file. Storing fewer or shorter reads is
the only thing that would shrink an archive substantially.

In memory the table is the same either way: `self.data` is a list of rows in
`self.fields` order, so nothing downstream of `read()` needs to know which
encoding a file used.
'''

import json

from .otu_table_entry import OtuTableEntry


class ArchiveOtuTable:
    version = 5

    FIELDS_VERSION1 = str.split(
        'gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known?'
    )
    FIELDS_VERSION2 = str.split(
        'gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences'
    )
    FIELDS_VERSION3 = str.split(
        'gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies'
    )
    FIELDS_VERSION4 = str.split(
        'gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method'
    )
    FIELDS_VERSION5 = FIELDS_VERSION4 + str.split('reads_with_repaired_deletions')
    FIELDS_OF_EACH_VERSION = [
        FIELDS_VERSION1,
        FIELDS_VERSION2,
        FIELDS_VERSION3,
        FIELDS_VERSION4,
        FIELDS_VERSION5,
    ]
    FIELDS = FIELDS_OF_EACH_VERSION[version - 1]

    # The first version written column-wise rather than row-wise.
    FIRST_COLUMNAR_VERSION = 5

    # The one field whose values are dereplicated into the shared 'reads' list,
    # being both by far the largest (~85% of an uncompressed archive) and the
    # only one that repeats across rows, since a read hitting two markers is
    # otherwise stored once per marker.
    UNALIGNED_SEQUENCE_FIELD_NAME = 'read_unaligned_sequences'

    # Reads are referred to by position in 'reads', but writing those positions
    # out as-is costs more compressed than the duplicate reads it saves: they are
    # mostly-ascending integers, which no compressor can do much with, whereas
    # the duplicate reads they replace were being matched almost for free. So a
    # read is instead referred to by NEXT_NEW_READ when it is the next one not
    # yet used, and by (position + 1) only when it is a repeat of an earlier one.
    # Since reads are numbered in the order the rows first use them, that turns
    # the common case into a run of identical values, which does compress away.
    NEXT_NEW_READ = 0

    READ_NAME_FIELD_INDEX = 6
    SAMPLE_ID_FIELD_INDEX = FIELDS_VERSION4.index('sample')
    UNALIGNED_SEQUENCE_FIELD_INDEX = FIELDS_VERSION4.index('read_unaligned_sequences')
    EQUAL_BEST_HIT_TAXONOMIES_INDEX = FIELDS_VERSION4.index('equal_best_hit_taxonomies')
    TAXONOMY_ASSIGNMENT_METHOD_INDEX = FIELDS_VERSION4.index('taxonomy_assignment_method')
    COVERAGE_FIELD_INDEX = FIELDS_VERSION4.index('coverage')
    TAXONOMY_FIELD_INDEX = FIELDS_VERSION4.index('taxonomy')
    NUCLEOTIDES_ALIGNED_FIELD_INDEX = FIELDS_VERSION4.index('nucleotides_aligned')
    TAXONOMY_BY_KNOWN_FIELD_INDEX = FIELDS_VERSION4.index('taxonomy_by_known?')
    READS_WITH_REPAIRED_DELETIONS_FIELD_INDEX = FIELDS_VERSION5.index('reads_with_repaired_deletions')

    def __init__(self, singlem_packages=None):
        self.singlem_packages = singlem_packages
        self.fields = self.FIELDS
        self.data = []
        if singlem_packages is not None:
            self.alignment_hmm_sha256s = list([s.alignment_hmm_sha256() for s in self.singlem_packages])
            self.singlem_package_sha256s = list([s.singlem_package_sha256() for s in self.singlem_packages])
        else:
            self.alignment_hmm_sha256s = None
            self.singlem_package_sha256s = None

    def add(self, new_otus):
        for otu in new_otus:
            self.data.append(otu.data)

    def write_to(self, output_io):
        json.dump(self._to_json(), output_io)

    def _to_json(self):
        j = {
            "version":
                self.version,
            "alignment_hmm_sha256s":
                self.alignment_hmm_sha256s
                if self.alignment_hmm_sha256s else [s.alignment_hmm_sha256() for s in self.singlem_packages],
            "singlem_package_sha256s":
                self.singlem_package_sha256s
                if self.singlem_package_sha256s else [s.singlem_package_sha256() for s in self.singlem_packages],
            'fields':
                self.fields,
        }
        if self.version >= ArchiveOtuTable.FIRST_COLUMNAR_VERSION:
            j.update(self._columnar_body())
        else:
            j["otus"] = self.data
        return j

    def _rows_padded_to_fields(self):
        '''self.data, with each row as long as self.fields.

        Rows can be shorter when entries parsed from an older archive are added
        to a table declaring a newer field list. The trailing fields are unknown
        rather than empty, so they become None.'''
        width = len(self.fields)
        rows = []
        for row in self.data:
            if len(row) > width:
                raise Exception("An OTU has %i fields but the table declares only %i" %
                                (len(row), width))
            rows.append(row if len(row) == width else list(row) + [None] * (width - len(row)))
        return rows

    @staticmethod
    def _is_constant_column(values):
        '''True if every value is the same scalar. Lists and dicts are never
        hoisted, being both unhashable and unlikely to repeat.'''
        first = values[0]
        if not isinstance(first, (str, int, float, bool)) and first is not None:
            return False
        # True == 1 in Python, so the type is compared too rather than letting a
        # bool column collapse onto an int one (or vice versa).
        return all(type(value) is type(first) and value == first for value in values)

    def _columnar_body(self):
        '''The version 5 encoding of self.data: columns, minus the fields that
        are constant across every row, with read sequences dereplicated.'''
        rows = self._rows_padded_to_fields()
        columns = {field: [row[i] for row in rows] for (i, field) in enumerate(self.fields)}

        constant_fields = {}
        if len(rows) > 0:
            for field in self.fields:
                if ArchiveOtuTable._is_constant_column(columns[field]):
                    constant_fields[field] = columns.pop(field)[0]

        reads = []
        if ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_NAME in columns:
            read_to_index = {}
            dereplicated = []
            for row_reads in columns[ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_NAME]:
                if row_reads is None:
                    dereplicated.append(None)
                    continue
                references = []
                for read in row_reads:
                    index = read_to_index.get(read)
                    if index is None:
                        read_to_index[read] = len(reads)
                        reads.append(read)
                        references.append(ArchiveOtuTable.NEXT_NEW_READ)
                    else:
                        references.append(index + 1)
                dereplicated.append(references)
            columns[ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_NAME] = dereplicated

        return {
            # Recorded explicitly because neither the columns nor the constants
            # give the row count when every field happens to be constant.
            "num_otus": len(rows),
            "constant_fields": constant_fields,
            "reads": reads,
            "otus": columns,
        }

    @staticmethod
    def _rows_from_columnar_body(j):
        fields = j['fields']
        num_otus = j['num_otus']
        constant_fields = j['constant_fields']
        columns = j['otus']
        reads = j['reads']

        if ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_NAME in columns:
            next_new_read = 0
            rehydrated = []
            for row_references in columns[ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_NAME]:
                if row_references is None:
                    rehydrated.append(None)
                    continue
                row_reads = []
                for reference in row_references:
                    if reference == ArchiveOtuTable.NEXT_NEW_READ:
                        index = next_new_read
                        next_new_read += 1
                    else:
                        index = reference - 1
                    if index < 0 or index >= len(reads):
                        raise Exception("Archive OTU table refers to a read that does not exist")
                    row_reads.append(reads[index])
                rehydrated.append(row_reads)
            columns[ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_NAME] = rehydrated

        rows = [[None] * len(fields) for _ in range(num_otus)]
        for (i, field) in enumerate(fields):
            if field in constant_fields:
                value = constant_fields[field]
                for row in rows:
                    row[i] = value
            else:
                if field not in columns:
                    raise Exception("Archive OTU table is missing the '%s' field" % field)
                column = columns[field]
                if len(column) != num_otus:
                    raise Exception("Archive OTU table field '%s' has %i values but %i OTUs" %
                                    (field, len(column), num_otus))
                for (row, value) in zip(rows, column):
                    row[i] = value
        return rows

    @staticmethod
    def read(input_io, min_version=None):
        otus = ArchiveOtuTable()
        j = json.load(input_io)
        if j['version'] not in range(1, ArchiveOtuTable.version + 1):
            raise Exception("Wrong OTU table version detected")
        otus.version = j['version']
        if min_version is not None and otus.version < min_version:
            raise InsufficientArchiveOtuTableVersionException("OTU table version is too old, required: %d, found: %d" %
                                                              (min_version, otus.version))

        otus.alignment_hmm_sha256s = j['alignment_hmm_sha256s']
        otus.singlem_package_sha256s = j['singlem_package_sha256s']

        otus.fields = j['fields']
        if otus.fields != ArchiveOtuTable.FIELDS_OF_EACH_VERSION[j['version'] - 1]:
            raise Exception("Unexpected archive OTU table format detected")

        if otus.version >= ArchiveOtuTable.FIRST_COLUMNAR_VERSION:
            otus.data = ArchiveOtuTable._rows_from_columnar_body(j)
        else:
            otus.data = j['otus']
        return otus

    def sort(self):
        self.data.sort()

    def __iter__(self):
        for d in self.data:
            e = ArchiveOtuTableEntry()
            e.marker = d[0]
            e.sample_name = d[1]
            e.sequence = d[2]
            e.count = d[3]
            e.coverage = d[4]
            e.taxonomy = d[5]
            e.data = d
            e.fields = self.fields
            yield e


class ArchiveOtuTableEntry(OtuTableEntry):

    def read_names(self):
        '''Return a list of read names for this OTU'''
        return self.data[ArchiveOtuTable.READ_NAME_FIELD_INDEX]

    def read_unaligned_sequences(self):
        return self.data[ArchiveOtuTable.UNALIGNED_SEQUENCE_FIELD_INDEX]

    def equal_best_hit_taxonomies(self):
        return self.data[ArchiveOtuTable.EQUAL_BEST_HIT_TAXONOMIES_INDEX]

    def taxonomy_assignment_method(self):
        return self.data[ArchiveOtuTable.TAXONOMY_ASSIGNMENT_METHOD_INDEX]
        
    def nucleotides_aligned(self):
        return self.data[ArchiveOtuTable.NUCLEOTIDES_ALIGNED_FIELD_INDEX]
    
    def taxonomy_by_known(self):
        return self.data[ArchiveOtuTable.TAXONOMY_BY_KNOWN_FIELD_INDEX]

    def reads_with_repaired_deletions(self):
        '''Indices into read_names() of the reads whose window carries an
        ambiguous base that frameshift repair inserted, as opposed to one that
        was already in the raw read. None when there are none, and when the OTU
        came from an archive written before version 5, which did not record
        it.'''
        if len(self.data) <= ArchiveOtuTable.READS_WITH_REPAIRED_DELETIONS_FIELD_INDEX:
            return None
        return self.data[ArchiveOtuTable.READS_WITH_REPAIRED_DELETIONS_FIELD_INDEX]
    
    def coverage(self):
        return self.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX]
    
    def taxonomy(self):
        return self.data[ArchiveOtuTable.TAXONOMY_FIELD_INDEX]


class InsufficientArchiveOtuTableVersionException(Exception):
    pass
