import os
import logging
import itertools
import shutil
import extern
import tempfile
import json
import pandas as pd

import zenodo_backpack
from zenodo_backpack import ZenodoBackpackMalformedException

from .singlem_package import SingleMPackage
from .sequence_classes import SeqReader
from .metapackage_read_name_store import MetapackageReadNameStore

# These constants should be updated in tandem
# When updating, also update the name in the set_env_vars.sh script in admin/
DATA_DEFAULT_VERSION = '5.4.0'
CURRENT_GTDB_DATABASE_VERSION = 'R226'

GTDB_DATABASE_NAME = 'Genome Taxonomy Database (GTDB)'
DATA_ENVIRONMENT_VARIABLE = 'SINGLEM_METAPACKAGE_PATH'
DATA_DOI = '10.5281/zenodo.5739611' # This is the series DOI, not the individual metapackage's DOI

CUSTOM_TAXONOMY_DATABASE_NAME = 'custom_taxonomy_database'

LYREBIRD_DATA_DEFAULT_VERSION = '0.3.1'
LYREBIRD_DATA_ENVIRONMENT_VARIABLE = 'LYREBIRD_METAPACKAGE_PATH'
LYREBIRD_DATA_DOI = '10.5281/zenodo.14768887'


class Metapackage:
    '''A class for a set of SingleM packages, plus prefilter DB'''

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    VERSION_KEY = 'singlem_metapackage_version'
    PREFILTER_DB_PATH_KEY = 'prefilter_db_path'
    SINGLEM_PACKAGES = 'singlem_packages'
    NUCLEOTIDE_SDB = 'nucleotide_sdb'
    SQLITE_DB_PATH_KEY = 'sqlite_db_path_key'
    TAXON_GENOME_LENGTHS_KEY = 'taxon_genome_lengths'
    TAXONOMY_DATABASE_NAME_KEY = 'taxonomy_database_name'
    TAXONOMY_DATABASE_VERSION_KEY = 'taxonomy_database_version'
    DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS_KEY = 'diamond_prefilter_performance_parameters'
    DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS = 'diamond_taxonomy_assignment_performance_parameters'
    MAKEIDX_SENSITIVITY_PARAMS = 'makeidx_sensitivity_params'
    AVG_NUM_GENES_PER_SPECIES = 'avg_num_genes_per_species'

    _CURRENT_FORMAT_VERSION = 6

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
                    '4': [
                        VERSION_KEY,
                        PREFILTER_DB_PATH_KEY,
                        NUCLEOTIDE_SDB,
                        SQLITE_DB_PATH_KEY,
                        TAXON_GENOME_LENGTHS_KEY,
                        ],
                    '5': [
                        VERSION_KEY,
                        PREFILTER_DB_PATH_KEY,
                        NUCLEOTIDE_SDB,
                        SQLITE_DB_PATH_KEY,
                        TAXON_GENOME_LENGTHS_KEY,
                        TAXONOMY_DATABASE_NAME_KEY,
                        TAXONOMY_DATABASE_VERSION_KEY,
                        DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS_KEY,
                        DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS,
                        MAKEIDX_SENSITIVITY_PARAMS,
                        ],
                    '6': [
                        VERSION_KEY,
                        PREFILTER_DB_PATH_KEY,
                        NUCLEOTIDE_SDB,
                        SQLITE_DB_PATH_KEY,
                        TAXON_GENOME_LENGTHS_KEY,
                        TAXONOMY_DATABASE_NAME_KEY,
                        TAXONOMY_DATABASE_VERSION_KEY,
                        DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS_KEY,
                        DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS,
                        MAKEIDX_SENSITIVITY_PARAMS,
                        AVG_NUM_GENES_PER_SPECIES,
                        ],
                      }

    def __init__(self, package_paths=None, prefilter_path=None):
        # Array of gpkg names to SingleMPackage objects
        self._hmms_and_positions = {}
        self._prefilter_path = prefilter_path

        if package_paths:
            self.singlem_packages = [SingleMPackage.acquire(path) for path in package_paths]
            logging.info("Loaded %i SingleM packages" % len(self.singlem_packages))
            for pkg in self.singlem_packages:
                self._hmms_and_positions[pkg.base_directory()] = pkg
        # Otherwise return an empty metapackage

    @staticmethod
    def acquire_default_backpack():
        logging.debug("Acquiring SingleM packages from environment variable")
        if DATA_ENVIRONMENT_VARIABLE not in os.environ:
            raise Exception("The {} environment variable, which points to the default data directory, is not set. To download the default SingleM metapackage, use 'singlem data'".format(DATA_ENVIRONMENT_VARIABLE))
        try:
            backpack = zenodo_backpack.acquire(env_var_name=DATA_ENVIRONMENT_VARIABLE, version=DATA_DEFAULT_VERSION)
        except KeyError:
            # Has the payload directory rather than the base directory been specified?
            original_directory = os.environ[DATA_ENVIRONMENT_VARIABLE]
            from pathlib import Path
            path = Path(original_directory)
            backpack = zenodo_backpack.acquire(path=path.parent, version=DATA_DEFAULT_VERSION)
        except ZenodoBackpackMalformedException as e:
            raise Exception("The metapackage defined by the {} environment variable is either malformed or does not match the version encoded in the version of SingleM installed ({}). If you are wanting to run a custom metapackage (or one newer than the installed software) then use the --metapackage flag, rather than specifying the metapackage through the {} environment variable.".format(
                DATA_ENVIRONMENT_VARIABLE,
                DATA_DEFAULT_VERSION,
                DATA_ENVIRONMENT_VARIABLE,
                )) from e
        return backpack

    @staticmethod
    def acquire_lyrebird_backpack():
        logging.debug("Acquiring Lyrebird packages from environment variable")
        if LYREBIRD_DATA_ENVIRONMENT_VARIABLE not in os.environ:
            raise Exception("The {} environment variable, which points to the default data directory, is not set. To download the default Lyrebird metapackage, use 'lyrebird data'".format(LYREBIRD_DATA_ENVIRONMENT_VARIABLE))
        try:
            backpack = zenodo_backpack.acquire(env_var_name=LYREBIRD_DATA_ENVIRONMENT_VARIABLE, version=LYREBIRD_DATA_DEFAULT_VERSION)
        except KeyError:
            # Has the payload directory rather than the base directory been specified?
            original_directory = os.environ[LYREBIRD_DATA_ENVIRONMENT_VARIABLE]
            from pathlib import Path
            path = Path(original_directory)
            backpack = zenodo_backpack.acquire(path=path.parent, version=LYREBIRD_DATA_DEFAULT_VERSION)
        except ZenodoBackpackMalformedException as e:
            raise Exception("The metapackage defined by the {} environment variable is either malformed or does not match the version encoded in the version of SingleM installed ({}). If you are wanting to run a custom metapackage (or one newer than the installed software) then use the --metapackage flag, rather than specifying the metapackage through the {} environment variable.".format(
                LYREBIRD_DATA_ENVIRONMENT_VARIABLE,
                LYREBIRD_DATA_DEFAULT_VERSION,
                LYREBIRD_DATA_ENVIRONMENT_VARIABLE,
                )) from e
        return backpack

    @staticmethod
    def acquire_default():
        '''Acquire the default metapackage'''
        backpack = Metapackage.acquire_default_backpack()
        return Metapackage.acquire(backpack.payload_directory_string())
    
    @staticmethod
    def acquire_lyrebird():
        '''Acquire the default Lyrebird metapackage'''
        backpack = Metapackage.acquire_lyrebird_backpack()
        return Metapackage.acquire(backpack.payload_directory_string())


    @staticmethod
    def acquire(metapackage_path):
        with open(os.path.join(
            metapackage_path,
            SingleMPackage._CONTENTS_FILE_NAME)) as f:

            contents_hash = json.load(f)

        if Metapackage.VERSION_KEY not in contents_hash:
            # If the user specifies a .zb directory, acquire the
            # payload_directory
            zb_version_key = 'zenodo_backpack_version'
            if zb_version_key in contents_hash:
                logging.info("Acquiring SingleM metapackage from Zenodo backpack directory specified ..")
                backpack = zenodo_backpack.acquire(
                    path = metapackage_path)
                return Metapackage.acquire(backpack.payload_directory_string())
            else:
                raise Exception("SingleM metapackage directory does not contain a {} key, and it also does not appear to be a ZenodoBackpack directory".format(Metapackage.VERSION_KEY))

        v=contents_hash[Metapackage.VERSION_KEY]
        logging.debug("Loading version %i SingleM metapackage: %s" % (v, metapackage_path))

        if v not in range(1,7):
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
            if contents_hash[Metapackage.NUCLEOTIDE_SDB] is not None:
                mpkg._nucleotide_sdb_path = os.path.join(metapackage_path, contents_hash[Metapackage.NUCLEOTIDE_SDB])
            else:
                # singlem metapackage was invoked with --no-nucleotide-sdb
                mpkg._nucleotide_sdb_path = None
        if v >= 4:
            if contents_hash[Metapackage.TAXON_GENOME_LENGTHS_KEY] is not None:
                mpkg._taxon_genome_lengths_path = os.path.join(metapackage_path, contents_hash[Metapackage.TAXON_GENOME_LENGTHS_KEY])
            else:
                # singlem metapackage was invoked with --no-taxon-genome-lengths
                mpkg._taxon_genome_lengths_path = None
        if v >= 5:
            mpkg._taxonomy_database_name = contents_hash[Metapackage.TAXONOMY_DATABASE_NAME_KEY]
            mpkg._taxonomy_database_version = contents_hash[Metapackage.TAXONOMY_DATABASE_VERSION_KEY]
            mpkg._diamond_prefilter_performance_parameters = contents_hash[Metapackage.DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS_KEY]
            mpkg._diamond_taxonomy_assignment_performance_parameters = contents_hash[Metapackage.DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS]
            mpkg._makeidx_sensitivity_params = contents_hash[Metapackage.MAKEIDX_SENSITIVITY_PARAMS]
        if v >= 6:
            if contents_hash.get(Metapackage.AVG_NUM_GENES_PER_SPECIES, None) is not None:
                mpkg._avg_num_genes_per_species = contents_hash.get(Metapackage.AVG_NUM_GENES_PER_SPECIES, None)
            else:
                mpkg._avg_num_genes_per_species = None
        return mpkg

    @staticmethod
    def download(**kwargs):
        '''Download a metapackage from Zenodo'''
        output_directory = kwargs.pop('output_directory', None)
        lyrebird = kwargs.pop('lyrebird', False)
        if not output_directory:
            raise Exception("Must specify an output directory to download a new default metapackage to")

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if lyrebird:
            data_version = LYREBIRD_DATA_DEFAULT_VERSION
            data_doi = LYREBIRD_DATA_DOI
            data_env_var = LYREBIRD_DATA_ENVIRONMENT_VARIABLE
        else:
            data_version = DATA_DEFAULT_VERSION
            data_doi = DATA_DOI
            data_env_var = DATA_ENVIRONMENT_VARIABLE
        
        logging.info(f"Downloading data version {data_version} with ZenodoBackpack from {data_doi} ..")
        backpack = zenodo_backpack.ZenodoBackpackDownloader().download_and_extract(
            output_directory,
            data_doi,
            progress_bar=True,
            # Unfortunately, we are stuck with zenodo_backpackge 0.2.0 because
            # the newer version cannot be deployed on conda-forge, due to this
            # open ticket. 
            # https://github.com/conda-forge/zenodo_backpack-feedstock/pull/1
            #
            # version=DATA_DEFAULT_VERSION 
            )
        logging.info("Finished downloading data")

        logging.info("The environment variable {} can now be set to {}".format(data_env_var, output_directory))
        logging.info("For instance, the following can be included in your .bashrc (requires logout and login after inclusion):")
        logging.info("export {}='{}'".format(data_env_var, os.path.abspath(backpack.base_directory)))

    @staticmethod
    def verify(**kwargs):
        '''Verify that the ZenodoBackpack is valid'''
        od = kwargs.pop('output_directory', None) # not used
        lyrebird = kwargs.pop('lyrebird', False)
        if od is not None:
            raise Exception("Verification of downloaded data does not require an output directory to be specified. Use the environment variable {} to specify the location of the downloaded data".format(DATA_ENVIRONMENT_VARIABLE))

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        if lyrebird:
            data_version = LYREBIRD_DATA_DEFAULT_VERSION
            logging.info("Acquiring Lyrebird packages from environment variable")
            if not LYREBIRD_DATA_ENVIRONMENT_VARIABLE in os.environ:
                raise Exception("The {} environment variable, which points to the default data directory, is not set. To download the default Lyrebird metapackage, use 'singlem data'".format(LYREBIRD_DATA_ENVIRONMENT_VARIABLE))
            backpack = Metapackage.acquire_lyrebird_backpack()
        else:
            data_version = DATA_DEFAULT_VERSION
            logging.info("Acquiring SingleM packages from environment variable")
            if not DATA_ENVIRONMENT_VARIABLE in os.environ:
                raise Exception("The {} environment variable, which points to the default data directory, is not set. To download the default SingleM metapackage, use 'singlem data'".format(DATA_ENVIRONMENT_VARIABLE))
            backpack = Metapackage.acquire_default_backpack()
        
        logging.info("Verifying data with ZenodoBackpack ..")
        zenodo_backpack.ZenodoBackpackDownloader().verify(backpack, passed_version=data_version)
        logging.info("Finished verifying data")


    @staticmethod
    def generate(**kwargs):
        singlem_packages = kwargs.pop('singlem_packages')
        nucleotide_sdb = kwargs.pop('nucleotide_sdb')
        prefilter_clustering_threshold = kwargs.pop('prefilter_clustering_threshold')
        output_path = kwargs.pop('output_path')
        threads = kwargs.pop('threads')
        prefilter_diamond_db = kwargs.pop('prefilter_diamond_db')
        taxon_genome_lengths_csv = kwargs.pop('taxon_genome_lengths')
        taxonomic_database_name = kwargs.pop('taxonomy_database_name')
        taxonomic_database_version = kwargs.pop('taxonomy_database_version')
        diamond_prefilter_performance_parameters = kwargs.pop('diamond_prefilter_performance_parameters')
        diamond_taxonomy_assignment_performance_parameters = kwargs.pop('diamond_taxonomy_assignment_performance_parameters')
        makeidx_sensitivity_params = kwargs.pop('makeidx_sensitivity_params')
        calculate_average_num_genes_per_species = kwargs.pop('calculate_average_num_genes_per_species', False)

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if calculate_average_num_genes_per_species not in (True, False):
            raise Exception("calculate_average_num_genes_per_species must be a boolean")

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
        
        taxonomy_marker_counts = None
        avg_num_genes_per_species = None

        # calculate average number of genes per species
        if calculate_average_num_genes_per_species:
            logging.debug("Calculating average number of genes per on-target species ..")
            total_count = 0
            taxonomy_marker_counts = {}
            for pkg in singlem_packages:
                spkg = SingleMPackage.acquire(pkg)
                tax_hash = spkg.taxonomy_hash()
                target_domains = spkg.target_domains()
                for seq_id, taxonomy in tax_hash.items():
                    if taxonomy[0].replace('d__','') in target_domains:
                        tax = 'Root;' + ';'.join(taxonomy) # add Root so condense doesn't break
                        if tax in taxonomy_marker_counts:
                            taxonomy_marker_counts[tax] += 1
                        else:
                            taxonomy_marker_counts[tax] = 1
                        total_count += 1
            avg_num_genes_per_species = total_count / len(taxonomy_marker_counts)
            logging.info("Average number of genes per species: {} for {} species".format(avg_num_genes_per_species, len(taxonomy_marker_counts)))

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

        # Copy taxon genome lengths csv into output directory
        if taxon_genome_lengths_csv:
            taxon_genome_lengths_csv_abspath = os.path.abspath(taxon_genome_lengths_csv)
            taxon_genome_lengths_csv_name = os.path.basename(taxon_genome_lengths_csv_abspath)
            taxon_genome_lengths_csv_path = os.path.join(output_path, taxon_genome_lengths_csv_name)
            logging.info("Copying taxon genome lengths csv {} to {} ..".format(taxon_genome_lengths_csv, taxon_genome_lengths_csv_path))
            shutil.copy(taxon_genome_lengths_csv, taxon_genome_lengths_csv_path)
        else:
            logging.info("Skipping taxon genome lengths csv")
            taxon_genome_lengths_csv_name = None

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
        if makeidx_sensitivity_params != None:
            diamond_prefilter_performance_parameters = diamond_prefilter_performance_parameters + " " + " ".join(makeidx_sensitivity_params.split())
            logging.info("Using DIAMOND makeidx parameters: {}".format(makeidx_sensitivity_params))
            extern.run('diamond makeidx --db {} {}'.format(prefilter_dmnd_path, makeidx_sensitivity_params))
        else:
            extern.run('diamond makeidx --db {}'.format(prefilter_dmnd_path))

        if not prefilter_diamond_db:
            os.remove(prefilter_path)

        logging.info("Generating read name taxonomy store ..")
        sqlitedb_path = os.path.join(output_path, 'read_taxonomies.sqlite3')
        MetapackageReadNameStore.generate(
            singlem_packages, sqlitedb_path, taxonomy_marker_counts)

        contents_hash = {Metapackage.VERSION_KEY: 6,
                        Metapackage.SINGLEM_PACKAGES: singlem_package_relpaths,
                        Metapackage.PREFILTER_DB_PATH_KEY: prefilter_dmnd_name,
                        Metapackage.NUCLEOTIDE_SDB: nucleotide_sdb_name,
                        Metapackage.SQLITE_DB_PATH_KEY: os.path.basename(sqlitedb_path),
                        Metapackage.TAXON_GENOME_LENGTHS_KEY: taxon_genome_lengths_csv_name,
                        Metapackage.TAXONOMY_DATABASE_NAME_KEY: taxonomic_database_name,
                        Metapackage.TAXONOMY_DATABASE_VERSION_KEY: taxonomic_database_version,
                        Metapackage.DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS_KEY: diamond_prefilter_performance_parameters,
                        Metapackage.DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS: diamond_taxonomy_assignment_performance_parameters,
                        Metapackage.MAKEIDX_SENSITIVITY_PARAMS: makeidx_sensitivity_params,
                        Metapackage.AVG_NUM_GENES_PER_SPECIES: avg_num_genes_per_species,
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
                        if "X" in seq:
                            continue
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

    def get_all_taxonomy_strings(self):
        store = MetapackageReadNameStore.acquire(self._sqlite_db_path)
        return store.get_all_taxonomy_strings()
    
    def get_taxon_marker_counts(self, taxons):
        store = MetapackageReadNameStore.acquire(self._sqlite_db_path)
        return store.get_marker_counts_of_species(taxons)
    
    def get_all_taxon_marker_counts(self):
        store = MetapackageReadNameStore.acquire(self._sqlite_db_path)
        return store.get_all_marker_counts()

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

    def taxon_genome_lengths(self):
        try:
            tsv = self._taxon_genome_lengths_path
        except AttributeError:
            # Happens when version < 3 or metapackage created from spkgs directly
            return None
        if tsv is None:
            return None
        return pd.read_csv(tsv, sep='\t')
    
    def taxonomy_database_name(self):
        try:
            return self._taxonomy_database_name
        except AttributeError:
            # Happens when version < 5 or metapackage created from spkgs directly
            return None

    def taxonomy_database_version(self):
        try:
            return self._taxonomy_database_version
        except AttributeError:
            # Happens when version < 5 or metapackage created from spkgs directly
            return None
    
    def diamond_prefilter_performance_parameters(self):
        try:
            return self._diamond_prefilter_performance_parameters
        except AttributeError:
            # Happens when version < 5 or metapackage created from spkgs directly
            return None

    def diamond_taxonomy_assignment_performance_parameters(self):
        try:
            return self._diamond_taxonomy_assignment_performance_parameters
        except AttributeError:
            # Happens when version < 5 or metapackage created from spkgs directly
            return None
    
    def makeidx_sensitivity_params(self):
        try:
            return self._makeidx_sensitivity_params
        except AttributeError:
            # Happens when version < 5 or metapackage created from spkgs directly
            return None
    
    def avg_num_genes_per_species(self):
        try:
            return self._avg_num_genes_per_species
        except AttributeError:
            # Happens when version < 6 or metapackage created from spkgs directly
            return None