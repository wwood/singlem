import re
import os
import csv
import logging
import itertools
import pkg_resources
import extern
import tempfile

from .singlem_package import SingleMPackage
from .sequence_classes import SeqReader

class Metapackage:
    '''A class for a set of SingleM packages'''

    def __init__(self, package_paths=None):
        # Array of gpkg names to SingleMPackage objects
        self._hmms_and_positions = {}

        if package_paths:
            self.singlem_packages = [SingleMPackage.acquire(path) for path in package_paths]
            logging.info("Loaded %i SingleM packages" % len(self.singlem_packages))
        else:
            # Prefer production DB directory
            pkg_resources_db_directory = 'data'

            pkg_paths = pkg_resources.resource_listdir('singlem',pkg_resources_db_directory)
            basedir = pkg_resources.resource_filename('singlem',pkg_resources_db_directory)
            logging.debug("Searching for SingleM packages via pkg_resources in %s .." % basedir)
            pkg_paths = [os.path.join(basedir,d) for d in pkg_paths if d[-5:]=='.spkg']
            if len(pkg_paths) == 0:
                raise Exception("Unable to find any SingleM packages using pkg_resources")

            logging.debug("Found %i SingleM packages: %s" % (len(pkg_paths),
                                                        ', '.join(pkg_paths)))
            self.singlem_packages = [SingleMPackage.acquire(path) for path in pkg_paths]

        for pkg in self.singlem_packages:
            self._hmms_and_positions[pkg.base_directory()] = pkg

    def get_dmnd(self):
        ''' Create temporary DIAMOND file for search method '''
        fasta_paths = [pkg.graftm_package().unaligned_sequence_database_path() for pkg in self.singlem_packages]
        temp_dmnd = tempfile.NamedTemporaryFile(mode="w", prefix='singlem-diamond-prefilter', 
                                                suffix='.dmnd', delete=False).name
        cmd = 'cat %s | '\
            'diamond makedb --in - --db %s' % (' '.join(fasta_paths), temp_dmnd)
        
        extern.run(cmd)
        
        return temp_dmnd
    
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

    def create_on_target_prefilter_fasta(self, output_prefilter_fasta_path):
        total_written_seqs = 0
        with open(output_prefilter_fasta_path,'w') as out:
            for pkg in self._hmms_and_positions.values():
                logging.info("Creating FASTA for {} ..".format(pkg.base_directory()))
                if pkg.version != 3:
                    raise Exception("Creating a prefilter DB only works on version 3 SingleM packages currently")
                tax_hash = pkg.graftm_package().taxonomy_hash()
                tax_hash_example_key = list(tax_hash.keys())[0]
                logging.debug("Read in {} taxonomies e.g. {}: {}".format(
                    len(tax_hash),
                    tax_hash_example_key,
                    tax_hash[tax_hash_example_key]
                ))
                seq_ids_to_write = set()
                total_seqs = 0
                for seq_id in pkg.get_sequence_ids():
                    total_seqs += 1
                    if tax_hash[seq_id][0].replace('d__','') in pkg.target_domains():
                        seq_ids_to_write.add(seq_id)
                logging.info("Found {} sequences in the target domain, out of {} total".format(
                    len(seq_ids_to_write), total_seqs))
                with open(pkg.graftm_package().unaligned_sequence_database_path()) as f:
                    for (name, seq, _) in SeqReader().readfq(f):
                        if name in seq_ids_to_write:
                            out.write(">{}\n{}\n".format(name, seq))
                            total_written_seqs += 1
        logging.info("Wrote {} sequences in total".format(total_written_seqs))
