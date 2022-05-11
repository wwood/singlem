import os
import logging
import itertools
import shutil
import extern
import tempfile
import json

import zenodo_backpack

from .singlem_package import SingleMPackage
from .sequence_classes import SeqReader
from .metapackage_read_name_store import MetapackageReadNameStore

DATA_DEFAULT_VERSION = '3.0.5'
DATA_ENVIRONMENT_VARIABLE = 'SINGLEM_METAPACKAGE_PATH'
DATA_DOI = '10.5281/zenodo.5739611'

class Metapackage:
    '''A class for a set of SingleM packages, plus prefilter DB'''

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    VERSION_KEY = 'singlem_metapackage_version'
    PREFILTER_DB_PATH_KEY = 'prefilter_db_path'
    SINGLEM_PACKAGES = 'singlem_packages'
    NUCLEOTIDE_SDB = 'nucleotide_sdb'
    SQLITE_DB_PATH_KEY = 'sqlite_db_path_key'

    _CURRENT_FORMAT_VERSION = 3

    _REQUIRED_KEYS = {'1': [
                            VERSION_KEY,
                            PREFILTER_DB_PATH_KEY,
                            ],
                    '3': [
                        VERSION_KEY,
                        PREFILTER_DB_PATH_KEY,
                        NUCLEOTIDE_SDB,
                        SQLITE_DB_PATH_KEY,
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
            raise Exception("No packages provided to create a new metapackage")
        for pkg in self.singlem_packages:
            self._hmms_and_positions[pkg.base_directory()] = pkg

    @staticmethod
    def acquire_default():
        '''Acquire the default metapackage'''
        logging.debug("Acquiring SingleM packages from environment variable")
        if not DATA_ENVIRONMENT_VARIABLE in os.environ:
            raise Exception("The {} environment variable, which points to the default data directory, is not set. To download the default SingleM metapackage, use 'singlem data'".format(DATA_ENVIRONMENT_VARIABLE))
        backpack = zenodo_backpack.acquire(env_var_name=DATA_ENVIRONMENT_VARIABLE, version=DATA_DEFAULT_VERSION)
        return Metapackage.acquire(backpack.payload_directory_string(enter_single_payload_directory=True))


    @staticmethod
    def acquire(metapackage_path):
        with open(os.path.join(
            metapackage_path,
            SingleMPackage._CONTENTS_FILE_NAME)) as f:

            contents_hash = json.load(f)

        v=contents_hash[Metapackage.VERSION_KEY]
        logging.debug("Loading version %i SingleM metapackage: %s" % (v, metapackage_path))

        if v not in (1,2,3):
            raise Exception("Bad SingleM metapackage version: %s" % str(v))

        spkg_relative_paths = contents_hash[Metapackage.SINGLEM_PACKAGES]

        mpkg = Metapackage(
            package_paths=[os.path.join(metapackage_path, pth) for pth in spkg_relative_paths],
            prefilter_path = os.path.join(metapackage_path, contents_hash[Metapackage.PREFILTER_DB_PATH_KEY]))
        mpkg._contents_hash = contents_hash
        mpkg._base_directory = metapackage_path
        mpkg.version = v

        if v >= 3:
            mpkg._sqlite_db_path = os.path.join(metapackage_path, contents_hash[Metapackage.SQLITE_DB_PATH_KEY])
            mpkg._nucleotide_sdb_path = os.path.join(metapackage_path, contents_hash[Metapackage.NUCLEOTIDE_SDB])

        return mpkg

    # @staticmethod
    # def acquire_from_common_arguments(**kwargs):
    #     '''
    #     Acquire a metapackage from a set of common command line arguments and
    #     environment variables.
    #     '''
    #     metapackage_path = kwargs.pop('metapackage_path', None)
    #     singlem_package_paths = kwargs.pop('singlem_packages', None)

    #     if metapackage_path and singlem_package_paths and singlem_package_paths != []:
    #         raise Exception("Cannot specify both a metapackage and singlem_packages")
    #     elif metapackage_path:
    #         return Metapackage.acquire(metapackage_path)
    #     elif not singlem_package_paths or singlem_package_paths == []:
    #         # Return the default set
    #         return Metapackage()
    #     else:
    #         return Metapackage(singlem_package_paths)

    @staticmethod
    def _grok_data_path(**kwargs):
        output_directory = kwargs.pop('output_directory', None)

        dir_from_env = os.environ.get(DATA_ENVIRONMENT_VARIABLE)

        if dir_from_env and output_directory:
            raise Exception("Cannot specify both an output directory and have a specified {}".format(DATA_ENVIRONMENT_VARIABLE))
        elif dir_from_env:
            output_directory = dir_from_env
        elif not output_directory:
            raise Exception("Either an output directory must be specified on the command line, or an environment variable must be set: {}".format(DATA_ENVIRONMENT_VARIABLE))
        
        return output_directory

    @staticmethod
    def download(**kwargs):
        '''Download a metapackage from Zenodo'''
        output_directory = Metapackage._grok_data_path(kwargs)
        kwargs.pop('output_directory', None)

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        logging.info("Downloading data with ZenodoBackpack ..")
        zenodo_backpack.ZenodoBackpackDownloader().download_and_extract(
            output_directory,
            DATA_DOI,
            progress_bar=True)
        logging.info("Finished downloading data")

        logging.info("The environment variable {} can now be set to {}".format(DATA_ENVIRONMENT_VARIABLE, output_directory))
        logging.info("For instance, the following can be included in your .bashrc (requires logout and login after inclusion):")
        logging.info("export {}='{}'".format(DATA_ENVIRONMENT_VARIABLE, output_directory))

    @staticmethod
    def verify(**kwargs):
        '''Verify that the ZenodoBackpack is valid'''
        output_path = Metapackage._grok_data_path(**kwargs)
        kwargs.pop('output_directory', None)

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        logging.info("Verifying data with ZenodoBackpack ..")
        zb = zenodo_backpack.acquire(output_path)
        zenodo_backpack.ZenodoBackpackDownloader().verify(zb, passed_version=DATA_DEFAULT_VERSION)
        logging.info("Finished verifying data")


    @staticmethod
    def generate(**kwargs):
        singlem_packages = kwargs.pop('singlem_packages')
        nucleotide_sdb = kwargs.pop('nucleotide_sdb')
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

        # Copy nucleotide SingleM db into output directory
        if nucleotide_sdb:
            # Remove trailing slash if present
            nucleotide_sdb_abspath = os.path.abspath(nucleotide_sdb)
            nucleotide_sdb_name = os.path.basename(nucleotide_sdb_abspath)
            nucleotide_sdb_path = os.path.join(output_path, nucleotide_sdb_name)
            logging.info("Copying SingleM db {} to {} ..".format(nucleotide_sdb, nucleotide_sdb_path))
            shutil.copytree(nucleotide_sdb, nucleotide_sdb_path)
        else:
            logging.info("Skipping SingleM db")
            nucleotide_sdb_name = None

        # Create on-target and dereplicated prefilter fasta file
        if prefilter_diamond_db:
            postfix = '.dmnd'
            if not prefilter_diamond_db.endswith(postfix):
                raise Exception("Predefined DIAMOND DB should end in ".format(postfix))
            prefilter_dmnd_name = os.path.basename(prefilter_diamond_db)
            prefilter_dmnd_path = os.path.join(output_path, prefilter_dmnd_name)
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

        logging.info("Generating read name taxonomy store ..")
        sqlitedb_path = os.path.join(output_path, 'read_taxonomies.sqlite3')
        MetapackageReadNameStore.generate(
            singlem_packages, sqlitedb_path)

        contents_hash = {Metapackage.VERSION_KEY: 3,
                        Metapackage.SINGLEM_PACKAGES: singlem_package_relpaths,
                        Metapackage.PREFILTER_DB_PATH_KEY: prefilter_dmnd_name,
                        Metapackage.NUCLEOTIDE_SDB: nucleotide_sdb_name,
                        Metapackage.SQLITE_DB_PATH_KEY: os.path.basename(sqlitedb_path),
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
                if pkg.version < 3:
                    raise Exception("Creating a prefilter DB only works on version 3 SingleM packages currently")
                tax_hash = pkg.taxonomy_hash()
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
        extern.run('cd-hit -n 3 -M 0 -i {} -o {} -T {} -c {}'.format(
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

    def describe(self):
        print('\t'.join(
            ('name','version','target_domains','gene_description')
        ))
        for spkg in self.singlem_packages:
            print('\t'.join((
                spkg.graftm_package_basename(),
                str(spkg.version),
                ','.join(spkg.target_domains()),
                spkg.gene_description()
            )))

    def get_taxonomy_of_reads(self, read_names):
        store = MetapackageReadNameStore.acquire(self._sqlite_db_path)
        return store.get_taxonomy_of_reads(read_names)

    def nucleotide_sdb(self):
        # import here so that we avoid tensorflow dependency if not needed
        from .sequence_database import SequenceDatabase
        db_path = self.nucleotide_sdb_path()
        if not db_path: return None
        return SequenceDatabase.acquire(db_path)
    
    def nucleotide_sdb_path(self):
        try:
            path = self._nucleotide_sdb_path
        except AttributeError:
            # Happens when version < 3 or metapackage created from spkgs directly
            return None
        return path
