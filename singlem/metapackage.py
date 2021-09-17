import re
import os
import csv
import logging
import itertools
import shutil
import pkg_resources
import extern
import tempfile
import json

from .singlem_package import SingleMPackage
from .sequence_classes import SeqReader

class Metapackage:
    '''A class for a set of SingleM packages, plus prefilter DB'''

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    VERSION_KEY = 'singlem_metapackage_version'
    PREFILTER_DB_PATH_KEY = 'prefilter_db_path'
    SINGLEM_PACKAGES = 'singlem_packages'

    _CURRENT_FORMAT_VERSION = 1

    _REQUIRED_KEYS = {'1': [
                            VERSION_KEY,
                            PREFILTER_DB_PATH_KEY,
                            ],
                      }

    def __init__(self, package_paths=None, prefilter_path=None):
        # Array of gpkg names to SingleMPackage objects
        self._hmms_and_positions = {}
        self._prefilter_path = prefilter_path

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

    @staticmethod
    def acquire(metapackage_path):
        with open(os.path.join(
            metapackage_path,
            SingleMPackage._CONTENTS_FILE_NAME)) as f:

            contents_hash = json.load(f)

        v=contents_hash[Metapackage.VERSION_KEY]
        logging.debug("Loading version %i SingleM metapackage: %s" % (v, metapackage_path))

        if v == 1:
            pkg = Metapackage()
        else:
            raise Exception("Bad SingleM metapackage version: %s" % str(v))

        spkg_relative_paths = contents_hash[Metapackage.SINGLEM_PACKAGES]

        mpkg = Metapackage(
            package_paths=[os.path.join(metapackage_path, pth) for pth in spkg_relative_paths],
            prefilter_path = os.path.join(metapackage_path, contents_hash[Metapackage.PREFILTER_DB_PATH_KEY]))
        mpkg._contents_hash = contents_hash
        mpkg._base_directory = metapackage_path
        return mpkg

    @staticmethod
    def generate(**kwargs):
        singlem_packages = kwargs.pop('singlem_packages')
        prefilter_clustering_threshold = kwargs.pop('prefilter_clustering_threshold')
        output_path = kwargs.pop('output_path')
        threads = kwargs.pop('threads')
        prefilter_diamond_db = kwargs.pop('prefilter_diamond_db')

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if os.path.exists(output_path):
            raise Exception("Not writing new SingleM metapackage to already existing file/directory with name %s" % output_path)
        os.mkdir(output_path)

        mpkg = Metapackage(package_paths=singlem_packages)

        # Copy singlem packages into output directory
        singlem_package_relpaths = []
        for pkg in singlem_packages:
            relpath = os.path.basename(os.path.abspath(pkg))
            if relpath in singlem_package_relpaths:
                raise Exception("Cannot have 2 singlem packages with the same basename in a metapackage. relpath was {}".format(relpath))
            singlem_package_relpaths.append(relpath)
            dest = os.path.join(output_path, relpath)
            logging.info("Copying package {} to be {} ..".format(pkg, dest))
            shutil.copytree(pkg, dest)

        # Create on-target and dereplicated prefilter fasta file
        if prefilter_diamond_db:
            if not prefilter_diamond_db.endswith('.dmnd'):
                raise Exception("Predefined DIAMOND DB should end in .dmnd")
            prefilter_dmnd_name = os.path.basename(prefilter_diamond_db)
            prefilter_dmnd_path = os.path.join(output_path, prefilter_name)
            shutil.copy(prefilter_diamond_db, prefilter_dmnd_path)
        else:
            prefilter_name = 'prefilter.fna'
            prefilter_path = os.path.join(output_path, prefilter_name)
            with tempfile.NamedTemporaryFile(prefix='singlem_metapackage') as t1:
                logging.info("Creating on-target FASTA file for prefilter ..")
                mpkg.create_on_target_prefilter_fasta(t1.name)

                logging.info("Dereplicating on target prefilter FASTA ..")
                mpkg.dereplicate_prefilter_fasta(
                    t1.name,
                    prefilter_path,
                    threads,
                    prefilter_clustering_threshold)

            # Create diamond DB indices of prefilter
            logging.info("Indexing prefilter DB with DIAMOND makedb ..")
            extern.run('diamond makedb --in {} --db {}.dmnd'.format(prefilter_path, prefilter_path))
            prefilter_dmnd_name = "{}.dmnd".format(prefilter_name)
            prefilter_dmnd_path = "{}.dmnd".format(prefilter_path)


        logging.info("Running DIAMOND makeidx of prefilter ..")
        extern.run('diamond makeidx --db {}'.format(prefilter_dmnd_path))

        if not prefilter_diamond_db:
            os.remove(prefilter_path)

        contents_hash = {Metapackage.VERSION_KEY: 1,
                        Metapackage.SINGLEM_PACKAGES: singlem_package_relpaths,
                        Metapackage.PREFILTER_DB_PATH_KEY: prefilter_dmnd_name
                        }

        # save contents file
        with open(os.path.join(output_path, Metapackage._CONTENTS_FILE_NAME), 'w') as f:
            json.dump(contents_hash, f)
        logging.info("Metapackage generated.")

    def prefilter_db_path(self):
        return self._prefilter_path

    def set_prefilter_db_path(self, path):
        self._prefilter_path = path

    def get_dmnd(self):
        ''' Create temporary DIAMOND file for search method '''
        fasta_paths = [pkg.graftm_package().unaligned_sequence_database_path() for pkg in self.singlem_packages]
        temp_dmnd = tempfile.NamedTemporaryFile(mode="w", prefix='singlem-diamond-prefilter',
                                                suffix='.dmnd', delete=False).name
        cmd = 'cat %s | '\
            'diamond makedb --in - --db %s' % (' '.join(fasta_paths), temp_dmnd)

        extern.run(cmd)
        extern.run("diamond makeidx -d {}".format(temp_dmnd))

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
                logging.info("Reading FASTA from {} ..".format(pkg.base_directory()))
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

    def dereplicate_prefilter_fasta(self, input_fasta_path, output_fasta_path, threads, clustering_threshold):
        '''Run CD-HIT to dereplicate the prefilter FASTA file'''
        # Use -n 3 since why not, and required for 0.6 threshold
        extern.run('cd-hit -n 3 -i {} -o {} -T {} -c {}'.format(
            input_fasta_path,
            output_fasta_path,
            threads,
            clustering_threshold
        ))
        os.remove("{}.clstr".format(output_fasta_path))
        count = 0
        with open(output_fasta_path) as f:
            for (name, seq, _) in SeqReader().readfq(f):
                count += 1
        logging.info("After clustering, {} sequences remained".format(count))
