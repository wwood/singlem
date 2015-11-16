import re
import os
import csv
import logging
from singlem_package import SingleMPackage

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
    def __init__(self):
        # Array of gpkg names to HmmAndPostion
        self.hmms_and_positions = {}
        db_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                       '..','db')
        pkg_paths = [d for d in os.listdir(db_directory) if d[-5:]=='.spkg']
        logging.debug("Found %i SingleM packages: %s", (len(pkg_paths),
                                                    ', '.join(pkg_paths)))
        if len(pkg_paths) == 0:
            raise Exception("Unable to find any SingleM packages in %s" % db_directory)
        self.singlem_packages = [SingleMPackage.acquire(os.path.join(db_directory, path)) for path in pkg_paths]

        for pkg in self.singlem_packages:
            self.hmms_and_positions[pkg.hmm_basename()] = \
                HmmAndPostion(pkg.graftm_package_path(),
                               pkg.hmm_path(),
                               pkg.singlem_position())

    def hmm_paths(self):
        'return an array of absolute paths to the hmms in this database'
        return [hp.hmm_filename for hp in self.hmms_and_positions.values()]

    def gpkg_basenames(self):
        return self.hmms_and_positions.keys()

    def gpkg_paths(self):
        return [h.gpkg_path for _, h in self.hmms_and_positions.iteritems()]

    def __iter__(self):
        for hp in self.hmms_and_positions.values():
            yield hp

class HmmAndPostion:
    def __init__(self, gpkg_path, hmm_filename, best_position):
        self.gpkg_path = gpkg_path
        self.hmm_filename = hmm_filename
        self.best_position = best_position

    def hmm_path(self):
        return os.path.join(self.gpkg_path, self.hmm_filename)
    
    def gpkg_basename(self):
        return os.path.basename(self.gpkg_path)
    
    def hmm_basename(self):
        return os.path.basename(self.hmm_filename)



