import json

from .otu_table_entry import OtuTableEntry


class ArchiveOtuTable:
    version = 4

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
    FIELDS_OF_EACH_VERSION = [
        FIELDS_VERSION1,
        FIELDS_VERSION2,
        FIELDS_VERSION3,
        FIELDS_VERSION4,
    ]
    FIELDS = FIELDS_OF_EACH_VERSION[version - 1]

    READ_NAME_FIELD_INDEX = 6
    SAMPLE_ID_FIELD_INDEX = FIELDS_VERSION4.index('sample')
    UNALIGNED_SEQUENCE_FIELD_INDEX = FIELDS_VERSION4.index('read_unaligned_sequences')
    EQUAL_BEST_HIT_TAXONOMIES_INDEX = FIELDS_VERSION4.index('equal_best_hit_taxonomies')
    TAXONOMY_ASSIGNMENT_METHOD_INDEX = FIELDS_VERSION4.index('taxonomy_assignment_method')
    COVERAGE_FIELD_INDEX = FIELDS_VERSION4.index('coverage')
    TAXONOMY_FIELD_INDEX = FIELDS_VERSION4.index('taxonomy')
    NUCLEOTIDES_ALIGNED_FIELD_INDEX = FIELDS_VERSION4.index('nucleotides_aligned')
    TAXONOMY_BY_KNOWN_FIELD_INDEX = FIELDS_VERSION4.index('taxonomy_by_known?')

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
        json.dump(
            {
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
                "otus":
                    self.data
            }, output_io)

    @staticmethod
    def read(input_io, min_version=None):
        otus = ArchiveOtuTable()
        j = json.load(input_io)
        if not j['version'] in [1, 2, 3, 4]:
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
    
    def coverage(self):
        return self.data[ArchiveOtuTable.COVERAGE_FIELD_INDEX]
    
    def taxonomy(self):
        return self.data[ArchiveOtuTable.TAXONOMY_FIELD_INDEX]


class InsufficientArchiveOtuTableVersionException(Exception):
    pass
