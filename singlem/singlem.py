import re
import os
import csv
import logging
from singlem_package import SingleMPackage
import itertools

class OrfMUtils:
    def un_orfm_name(self, name):
        return re.sub('_\d+_\d+_\d+$', '', name)

class TaxonomyFile:
    def __init__(self, taxonomy_file_path):
        self.sequence_to_taxonomy = {}
        utils = OrfMUtils()
        with open(taxonomy_file_path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.sequence_to_taxonomy[\
                      utils.un_orfm_name(row[0])] = row[1]

    def __getitem__(self, item):
        return self.sequence_to_taxonomy[item]

class HmmDatabase:
    def __init__(self, package_paths=None):
        # Array of gpkg names to SingleMPackage objects
        self._hmms_and_positions = {}

        if package_paths:
            self.singlem_packages = [SingleMPackage.acquire(path) for path in package_paths]
            logging.info("Loaded %i SingleM packages" % len(self.singlem_packages))
        else:
            # Prefer production DB directory
            db_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'db')
            if not os.path.exists(db_directory):
                # If no production DB directory exists, use the dev DB directory
                db_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                            '..','db')
            pkg_paths = [d for d in os.listdir(db_directory) if d[-5:]=='.spkg']
            logging.debug("Found %i SingleM packages: %s" % (len(pkg_paths),
                                                        ', '.join(pkg_paths)))
            if len(pkg_paths) == 0:
                raise Exception("Unable to find any SingleM packages in %s" % db_directory)
            self.singlem_packages = [SingleMPackage.acquire(os.path.join(db_directory, path)) for path in pkg_paths]

        for pkg in self.singlem_packages:
            self._hmms_and_positions[pkg.base_directory()] = pkg

    def protein_packages(self):
        return [pkg for pkg in self._hmms_and_positions.values() if pkg.is_protein_package()]

    def nucleotide_packages(self):
        return [pkg for pkg in self._hmms_and_positions.values() if not pkg.is_protein_package()]

    def protein_search_hmm_paths(self):
        'return an array of absolute paths to the protein hmms in this database'
        return list(itertools.chain(
            *[pkg.graftm_package().search_hmm_paths() for pkg in self.protein_packages()]))

    def nucleotide_search_hmm_paths(self):
        'return an array of absolute paths to the protein hmms in this database'
        return list(itertools.chain(
            *[pkg.graftm_package().search_hmm_paths() for pkg in self.nucleotide_packages()]))

    def __iter__(self):
        for hp in self._hmms_and_positions.values():
            yield hp

