from string import split
import json
class ArchiveOtuTable:
    version = 1
    def __init__(self, singlem_packages):
        self.singlem_packages = singlem_packages
        self.fields = split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned')
        self.data = []
        
    def write_to(self, output_io):
        json.dump({"version": self.version,
             "alignment_hmm_sha256s": [s.alignment_hmm_sha256() for s in self.singlem_packages],
             "singlem_package_sha256s": [s.singlem_package_sha256() for s in self.singlem_packages],
             'fields': self.fields,
             "otus": self.data},
                  output_io)
