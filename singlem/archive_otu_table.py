import json

from .otu_table_entry import OtuTableEntry


class ArchiveOtuTable:
    version = 3

    FIELDS_VERSION1 = str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known?')
    FIELDS_VERSION2 = str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences')
    FIELDS_VERSION3 = str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies')
    FIELDS = FIELDS_VERSION2

    READ_NAME_FIELD_INDEX=6

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

    def write_to(self, output_io):
        json.dump({"version": self.version,
             "alignment_hmm_sha256s": self.alignment_hmm_sha256s if self.alignment_hmm_sha256s else [s.alignment_hmm_sha256() for s in self.singlem_packages],
             "singlem_package_sha256s": self.singlem_package_sha256s if self.singlem_package_sha256s else [s.singlem_package_sha256() for s in self.singlem_packages],
             'fields': self.fields,
             "otus": self.data},
                  output_io)

    @staticmethod
    def read(input_io):
        otus = ArchiveOtuTable()
        j = json.load(input_io)
        if not j['version'] in [1,2,3]:
            raise Exception("Wrong OTU table version detected")
        otus.version = j['version']

        otus.alignment_hmm_sha256s = j['alignment_hmm_sha256s']
        otus.singlem_package_sha256s = j['singlem_package_sha256s']

        otus.fields = j['fields']
        if otus.fields != [ArchiveOtuTable.FIELDS_VERSION1,ArchiveOtuTable.FIELDS_VERSION2][j['version']-1]:
            raise Exception("Unexpected archive OTU table format detected")

        otus.data = j['otus']
        return otus

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
