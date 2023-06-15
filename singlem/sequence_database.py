import os
import tempfile
import logging
import subprocess
import sqlite3
import glob
import json
import itertools
import csv
import extern
import numpy as np

from sqlalchemy import create_engine, select, distinct

import Bio.Data.CodonTable

from .otu_table import OtuTable
from .singlem_database_models import NucleotideSequence, NucleotidesProteins, ProteinSequence, Taxonomy, Otu, Marker
from .otu_table_entry import OtuTableEntry

DEFAULT_NUM_THREADS = 1

ANNOY_INDEX_FORMAT = 'annoy'
NMSLIB_INDEX_FORMAT = 'nmslib'
SCANN_INDEX_FORMAT = 'scann'
SCANN_NAIVE_INDEX_FORMAT = 'scann-naive'
SMAFA_NAIVE_INDEX_FORMAT = 'smafa-naive'

NUCLEOTIDE_DATABASE_TYPE = 'nucleotide'
PROTEIN_DATABASE_TYPE = 'protein'

class SequenceDatabase:
    version = 5
    SQLITE_DB_NAME = 'otus.sqlite3'
    _marker_to_nmslib_nucleotide_index_file = {}
    _marker_to_nmslib_protein_index_file = {}
    _marker_to_annoy_nucleotide_index_file = {}
    _marker_to_annoy_protein_index_file = {}
    _marker_to_scann_nucleotide_index_file = {}
    _marker_to_scann_protein_index_file = {}
    _marker_to_scann_naive_nucleotide_index_file = {}
    _marker_to_scann_naive_protein_index_file = {}
    _marker_to_smafa_naive_nucleotide_index_file = {}

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    VERSION_KEY = 'singlem_database_version'

    _REQUIRED_KEYS = {5: [VERSION_KEY]}

    NUCLEOTIDE_TYPE = 'nucleotide'
    PROTEIN_TYPE = 'protein'

    _marker_cache = None
    _taxonomy_cache = None

    def get_marker_via_cache(self, marker_id):
        if self._marker_cache is None:
            logging.info('Loading marker cache')
            with self.engine.connect() as conn:
                self._marker_cache = Marker.generate_python_index(conn)
        return self._marker_cache[marker_id]

    def get_taxonomy_via_cache(self, taxonomy_id):
        if self._taxonomy_cache is None:
            logging.info('Loading taxonomy cache')
            with self.engine.connect() as conn:
                self._taxonomy_cache = Taxonomy.generate_python_index(conn)
        return self._taxonomy_cache[taxonomy_id]

    def add_sequence_db(self, marker_name, db_path, index_format, sequence_type):
        if index_format == NMSLIB_INDEX_FORMAT:
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                self._marker_to_nmslib_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                self._marker_to_nmslib_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format == ANNOY_INDEX_FORMAT:
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                self._marker_to_annoy_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                self._marker_to_annoy_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format == SCANN_INDEX_FORMAT:
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                self._marker_to_scann_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                self._marker_to_scann_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format == SCANN_NAIVE_INDEX_FORMAT:
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                self._marker_to_scann_naive_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                self._marker_to_scann_naive_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format == SMAFA_NAIVE_INDEX_FORMAT:
            
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                self._marker_to_smafa_naive_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
                raise NotImplementedError()
                # self._marker_to_smafa_naive_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        else:
            raise Exception("Unknown index type {}".format(index_format))

    def get_sequence_index(self, marker_name, index_format, sequence_type):
        index = None
        if index_format == NMSLIB_INDEX_FORMAT:
            if sequence_type == 'nucleotide':
                if marker_name in self._marker_to_nmslib_nucleotide_index_file:
                    index_path = self._marker_to_nmslib_nucleotide_index_file[marker_name]
                    index = SequenceDatabase._nucleotide_nmslib_init()
                    logging.debug("Loading index for {} from {}".format(marker_name, index_path))
                    index.loadIndex(index_path, load_data=True)
                    return index
            elif sequence_type == 'protein':
                if marker_name in self._marker_to_nmslib_protein_index_file:
                    index_path = self._marker_to_nmslib_protein_index_file[marker_name]
                    index = SequenceDatabase._nucleotide_nmslib_init()
                    logging.debug("Loading index for {} from {}".format(marker_name, index_path))
                    index.loadIndex(index_path, load_data=True)
                    return index
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format == ANNOY_INDEX_FORMAT:
            if sequence_type == 'nucleotide':
                if marker_name in self._marker_to_annoy_nucleotide_index_file:
                    index_path = self._marker_to_annoy_nucleotide_index_file[marker_name]
                    index = self._nucleotide_annoy_init()
                    logging.debug("Loading index for {} from {}".format(marker_name, index_path))
                    index.load(index_path)
                    return index
            elif sequence_type == 'protein':
                if marker_name in self._marker_to_annoy_protein_index_file:
                    index_path = self._marker_to_annoy_protein_index_file[marker_name]
                    index = self._protein_annoy_init()
                    logging.debug("Loading index for {} from {}".format(marker_name, index_path))
                    index.load(index_path)
                    return index
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format in [SCANN_INDEX_FORMAT, SCANN_NAIVE_INDEX_FORMAT]:
            import tensorflow as tf
            if sequence_type == 'nucleotide':
                if index_format == SCANN_INDEX_FORMAT:
                    markers_to_paths = self._marker_to_scann_nucleotide_index_file
                else:
                    markers_to_paths = self._marker_to_scann_naive_nucleotide_index_file
            elif sequence_type == 'protein':
                if index_format == SCANN_INDEX_FORMAT:
                    markers_to_paths = self._marker_to_scann_protein_index_file
                else:
                    markers_to_paths = self._marker_to_scann_naive_protein_index_file
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)

            import scann # only load when needed to speed start-up
            if marker_name in markers_to_paths:
                index_path = markers_to_paths[marker_name]
                logging.debug("Loading index for {} from {}".format(marker_name, index_path))
                if index_format == SCANN_INDEX_FORMAT:
                    # Can fail when there are too few sequences?
                    index = scann.scann_ops_pybind.load_searcher(index_path)
                else:
                    reloaded = tf.compat.v2.saved_model.load(export_dir=index_path)
                    index = scann.scann_ops.searcher_from_module(reloaded)
                return index
        elif index_format == SMAFA_NAIVE_INDEX_FORMAT:
            if sequence_type == 'nucleotide':
                if marker_name in self._marker_to_smafa_naive_nucleotide_index_file:
                    index_path = self._marker_to_smafa_naive_nucleotide_index_file[marker_name]
                    return index_path
            elif sequence_type == 'protein':
                raise NotImplementedError("SMAFA naive protein index not implemented")
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        else:
            raise Exception("Unknown index type {}".format(index_format))
        
        logging.warn("No %s / %s sequence index DB found for %s" % (index_format, sequence_type, marker_name))
        return None

    @staticmethod
    def _nucleotide_nmslib_init():
        import nmslib  # optional dependency
        return nmslib.init(space='bit_hamming', data_type=nmslib.DataType.OBJECT_AS_STRING, dtype=nmslib.DistType.INT, method='hnsw')

    @staticmethod
    def _protein_nmslib_init():
        import nmslib  # optional dependency
        return nmslib.init(space='bit_hamming', data_type=nmslib.DataType.OBJECT_AS_STRING, dtype=nmslib.DistType.INT, method='hnsw')

    def _nucleotide_annoy_init(self):
        example_seq = self.sqlalchemy_connection.execute(select(NucleotideSequence).limit(1)).first().sequence
        ndim = len(example_seq)*5
        from annoy import AnnoyIndex
        return AnnoyIndex(ndim, 'hamming')

    def _protein_annoy_init(self):
        example_seq = self.sqlalchemy_connection.execute(select(ProteinSequence).limit(1)).first().protein_sequence
        ndim = len(example_seq)*len(AA_ORDER)
        from annoy import AnnoyIndex
        return AnnoyIndex(ndim, 'hamming')

    @staticmethod
    def acquire(path, min_version=None):
        db = SequenceDatabase()
        db.base_directory = path

        contents_path = os.path.join(
            path, SequenceDatabase._CONTENTS_FILE_NAME)
        if not os.path.exists(contents_path):
            logging.error("The SingleM database located at '{}' did not contain a contents file.".format(
                path) +
                          "This means that the DB is not at that location, is corrupt, or was generated by SingleM version 0.11.0 or older. So unfortunately the DB could not be loaded.")
            raise Exception("Failed to find contents file in SingleM DB {}".format(
                contents_path))
        with open(contents_path) as f:
            db._contents_hash = json.load(f)

        found_version = db._contents_hash[SequenceDatabase.VERSION_KEY]
        logging.debug("Loading version {} SingleM database: {}".format(
            found_version, path))
        if min_version is not None:
            if found_version < min_version:
                raise Exception("SingleM database version {} is too old for this operations. Version {} or newer is required.".format(
                    found_version, min_version))
        if found_version == 5:
            for key in SequenceDatabase._REQUIRED_KEYS[found_version]:
                if key not in db._contents_hash:
                    raise Exception(
                        "Unexpectedly did not find required key {} in SingleM database contents file: {}".format(
                            key, path))
        else:
            raise Exception("Unexpected SingleM DB version found: {}".format(found_version))

        db.sqlite_file = os.path.join(path, SequenceDatabase.SQLITE_DB_NAME)
        db.engine = create_engine("sqlite:///" + db.sqlite_file)
        db.sqlalchemy_connection = db.engine.connect()

        nmslib_nucleotide_index_files = glob.glob("%s/nucleotide_indices_nmslib/*.nmslib_index" % path)
        logging.debug("Found nmslib_nucleotide_index_files: %s" % ", ".join(nmslib_nucleotide_index_files))
        for g in nmslib_nucleotide_index_files:
            marker = os.path.basename(g).replace('.nmslib_index','')
            db.add_sequence_db(marker, g, NMSLIB_INDEX_FORMAT,'nucleotide')

        nmslib_protein_index_files = glob.glob("%s/protein_indices_nmslib/*.nmslib_index" % path)
        logging.debug("Found nmslib_protein_index_files: %s" % ", ".join(nmslib_protein_index_files))
        for g in nmslib_protein_index_files:
            marker = os.path.basename(g).replace('.nmslib_index','')
            db.add_sequence_db(marker, g, NMSLIB_INDEX_FORMAT,'protein')

        annoy_nucleotide_index_files = glob.glob("%s/nucleotide_indices_annoy/*.annoy_index" % path)
        logging.debug("Found annoy_nucleotide_index_files: %s" % ", ".join(annoy_nucleotide_index_files))
        for g in annoy_nucleotide_index_files:
            marker = os.path.basename(g).replace('.annoy_index','')
            db.add_sequence_db(marker, g, ANNOY_INDEX_FORMAT,'nucleotide')

        annoy_protein_index_files = glob.glob("%s/protein_indices_annoy/*.annoy_index" % path)
        logging.debug("Found annoy_protein_index_files: %s" % ", ".join(nmslib_protein_index_files))
        for g in annoy_protein_index_files:
            marker = os.path.basename(g).replace('.annoy_index','')
            db.add_sequence_db(marker, g, ANNOY_INDEX_FORMAT,'protein')
        
        scann_nucleotide_index_files = glob.glob("%s/nucleotide_indices_scann/*" % path)
        logging.debug("Found scann_nucleotide_index_files: %s" % ", ".join(scann_nucleotide_index_files))
        for g in scann_nucleotide_index_files:
            marker = os.path.basename(g)
            db.add_sequence_db(marker, g, SCANN_INDEX_FORMAT,'nucleotide')
        
        scann_naive_nucleotide_index_files = glob.glob("%s/nucleotide_indices_scann_brute_force/*" % path)
        logging.debug("Found scann_brute_force_nucleotide_index_files: %s" % ", ".join(scann_naive_nucleotide_index_files))
        for g in scann_naive_nucleotide_index_files:
            marker = os.path.basename(g)
            db.add_sequence_db(marker, g, SCANN_NAIVE_INDEX_FORMAT,'nucleotide')
        
        scann_protein_index_files = glob.glob("%s/protein_indices_scann/*" % path)
        logging.debug("Found scann_protein_index_files: %s" % ", ".join(scann_protein_index_files))
        for g in scann_protein_index_files:
            marker = os.path.basename(g)
            db.add_sequence_db(marker, g, SCANN_INDEX_FORMAT,'protein')
        
        scann_naive_protein_index_files = glob.glob("%s/protein_indices_scann_brute_force/*" % path)
        logging.debug("Found scann_brute_force_protein_index_files: %s" % ", ".join(scann_naive_protein_index_files))
        for g in scann_naive_protein_index_files:
            marker = os.path.basename(g)
            db.add_sequence_db(marker, g, SCANN_NAIVE_INDEX_FORMAT,'protein')

        smafa_naive_nucleotide_index_files = glob.glob("%s/nucleotide_indices_smafa_naive/*" % path)
        logging.debug("Found smafa-naive: %s" % ", ".join(smafa_naive_nucleotide_index_files))
        for g in smafa_naive_nucleotide_index_files:
            marker = os.path.basename(g.replace('.smafa_naive_index',''))
            logging.debug("Found marker: %s" % marker)
            db.add_sequence_db(marker, g, SMAFA_NAIVE_INDEX_FORMAT,'nucleotide')

        return db

    @staticmethod
    def _grouper(iterable, n):
        args = [iter(iterable)] * n
        return itertools.zip_longest(*args, fillvalue=None)

    @staticmethod
    def create_from_otu_table(
        db_path,
        otu_table_collection, 
        num_threads=DEFAULT_NUM_THREADS,
        pregenerated_sqlite3_db=None,
        tmpdir=None,
        num_annoy_nucleotide_trees = 10, # ntrees are currently guesses
        num_annoy_protein_trees = 10,
        sequence_database_methods = [SMAFA_NAIVE_INDEX_FORMAT],
        sequence_database_types = [NUCLEOTIDE_DATABASE_TYPE]):

        if num_threads is None:
            num_threads = DEFAULT_NUM_THREADS

        if sequence_database_types is None:
            sequence_database_types = []
        for db_type in sequence_database_types:
            if db_type not in [NUCLEOTIDE_DATABASE_TYPE, PROTEIN_DATABASE_TYPE]:
                raise Exception("Unexpected sequence database type: {}".format(db_type))

        # ensure db does not already exist
        if os.path.exists(db_path):
            raise Exception("Cowardly refusing to overwrite already-existing database path '%s'" % db_path)
        logging.info("Creating SingleM database at {}".format(db_path))
        os.makedirs(db_path)

        # Create contents file
        contents_file_path = os.path.join(db_path, SequenceDatabase._CONTENTS_FILE_NAME)
        with open(contents_file_path, 'w') as f:
            json.dump({
                SequenceDatabase.VERSION_KEY: 5,
            }, f)

        if pregenerated_sqlite3_db:
            logging.info("Re-using previous SQLite database {}".format(pregenerated_sqlite3_db))
            sqlite_db_path = pregenerated_sqlite3_db

            # Symlink instead of cp to the old otus.sqlite3 because the main use
            # case for this is messing with the nmslib/annoy DB creation
            # parameters.
            from pathlib import Path
            Path(os.path.join(db_path,'otus.sqlite3')).symlink_to(os.path.abspath(pregenerated_sqlite3_db))

        else:
            # Dumping the table into SQL and then modifying it form there is
            # taking too long (or I don't understand SQL well enough). The main
            # issue is that creating a table with (id, sequence) take a while to
            # construct, and then a while to insert the foreign keys into the
            # main table. So instead take a more hacky approach and generate
            # them through GNU sort, which is multi-threaded.

            with tempfile.TemporaryDirectory() as my_tempdir:
                if tmpdir is not None:
                    my_tempdir = tmpdir

                total_otu_count = 0
                
                #################################################################
                # sort by marker and sequence, for later getting the unique
                # sequences, and collect the taxonomy hash.
                taxonomy_name_to_id = {}
                next_taxonomy_id = 1
                taxonomy_loading_path = os.path.join(my_tempdir, 'taxonomy.tsv')
                sorted_path = os.path.join(my_tempdir,'makedb_sort_output')
                # Have to set LC_COLLATE=C because otherwise dashes in sequences
                # can be ignored. See
                # https://serverfault.com/questions/95579/unix-sort-treats-dash-characters-as-invisible
                proc = subprocess.Popen(['bash','-c','LC_COLLATE=C sort --parallel={} --buffer-size=20% > {}'.format(num_threads, sorted_path)],
                    stdin=subprocess.PIPE,
                    stdout=None,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
                with open(taxonomy_loading_path, 'w') as taxonomy_loading_file:
                    for entry in otu_table_collection:
                        total_otu_count += 1

                        if entry.taxonomy in taxonomy_name_to_id:
                            taxonomy_id = taxonomy_name_to_id[entry.taxonomy]
                        else:
                            taxonomy_id = next_taxonomy_id
                            taxonomy_name_to_id[entry.taxonomy] = taxonomy_id
                            next_taxonomy_id += 1
                            taxonomy_loading_file.write('{}\t{}\n'.format(taxonomy_id, entry.taxonomy))

                        print("\t".join((entry.marker, entry.sequence, entry.sample_name, str(entry.count),
                                    str(entry.coverage), str(taxonomy_id))), file=proc.stdin)
                logging.info("Sorting {} OTU observations ..".format(total_otu_count))
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise Exception("Sort command returned non-zero exit status %i.\n"\
                        "STDERR was: %s" % (
                            proc.returncode, proc.stderr.read()))
                proc.stderr.close()

                logging.info("Loading taxonomy table ..")
                sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE taxonomy (id INTEGER PRIMARY KEY AUTOINCREMENT,"
                        " taxonomy text);\n"
                        '.separator "\\t"\n'
                        ".import {} taxonomy".format(taxonomy_loading_path))
                logging.info("Loaded {} taxonomy entries".format(len(taxonomy_name_to_id)))


                #################################################################
                logging.info("Creating numbered OTU table tsv")
                marker_index = 0
                sequence_index = 0
                numbered_table_file = os.path.join(my_tempdir,'makedb_numbered_output.tsv')
                # We ultimately sort the otus table by sample_name not sequence,
                # because that makes for faster dumping/extraction of all OTUs
                # from a sample. Use awk to add IDs to the OTU table after
                # sorting.
                cmd = "LC_COLLATE=C sort --parallel={} --buffer-size=20% |awk -F'\t' '{{OFS = FS}} {{print NR,$0}}' > {}".format(num_threads, numbered_table_file)
                proc = subprocess.Popen(['bash','-o','pipefail','-c',cmd],
                    stdin=subprocess.PIPE,
                    stdout=None,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
                numbered_marker_and_sequence_file = os.path.join(my_tempdir,'number_and_sequence_file')
                numbered_markers_file = os.path.join(my_tempdir,'makedb_numbered_markers_file')
                numbered_sequences_file = os.path.join(my_tempdir,'makedb_numbered_sequences_file')
                with open(sorted_path) as csvfile_in:
                    reader = csv.reader(csvfile_in, delimiter="\t")
                    last_marker = None
                    last_sequence = None
                    last_marker_wise_sequence_id = None
                    with open(numbered_markers_file,'w') as markers_output_table_io:
                        with open(numbered_sequences_file,'w') as nucleotides_output_table_io:
                            for row in reader:
                                new_nucleotide = False
                                if last_marker != row[0]:
                                    last_marker_wise_sequence_id = 0 # Annoy starts at 0, so go with that. nmslib is arbitrary.
                                    last_marker = row[0]
                                    marker_index += 1
                                    last_sequence = row[1]
                                    sequence_index += 1
                                    # marker id, marker name
                                    print("\t".join((str(marker_index), row[0])), file=markers_output_table_io)
                                    new_nucleotide = True
                                elif last_sequence != row[1]:
                                    last_marker_wise_sequence_id += 1
                                    last_sequence = row[1]
                                    sequence_index += 1
                                    new_nucleotide = True
                                print("\t".join(row[2:]+[str(marker_index),str(sequence_index),row[1], str(last_marker_wise_sequence_id)]), file=proc.stdin)
                                if new_nucleotide:
                                    # seuence_id, marker_id, sequence, marker_wise_sequnce_id
                                    print("\t".join((str(sequence_index), str(marker_index), row[1], str(last_marker_wise_sequence_id))), file=nucleotides_output_table_io)
                logging.info("Sorting OTU observations by sample_name ..")
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise Exception("Sort command returned non-zero exit status %i.\n"\
                        "STDERR was: %s" % (
                            proc.returncode, proc.stderr.read()))
                proc.stderr.close()

                #################################################################
                logging.info("Importing OTU table into SQLite ..")
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE otus (id INTEGER PRIMARY KEY AUTOINCREMENT,"
                        " sample_name text, num_hits int, coverage float, taxonomy_id int, marker_id int, sequence_id int, sequence text, marker_wise_sequence_id int);\n"
                        '.separator "\\t"\n'
                        ".import {} otus".format(numbered_table_file))

                logging.info("Importing markers table into SQLite ..")
                sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE markers (id INTEGER PRIMARY KEY,"
                        " marker text);\n"
                        '.separator "\\t"\n'
                        ".import {} markers".format(numbered_markers_file))

                logging.info("Importing sequence table into SQLite ..")
                sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE nucleotides (id INTEGER PRIMARY KEY,"
                        " marker_id int, sequence text, marker_wise_id int);\n"
                        '.separator "\\t"\n'
                        ".import {} nucleotides".format(numbered_sequences_file))


                logging.info("Creating SQL indexes on otus ..")
                sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
                logging.debug("Connecting to db %s" % sqlite_db_path)
                db = sqlite3.connect(sqlite_db_path)
                c = db.cursor()

                c.execute("CREATE INDEX otu_sample_name on otus (sample_name)")
                c.execute("CREATE INDEX otu_taxonomy on otus (taxonomy_id)")
                c.execute("CREATE INDEX otu_marker on otus (marker_id)")
                c.execute("CREATE INDEX otu_sequence on otus (sequence_id)")
                c.execute("CREATE INDEX otu_marker_wise_sequence_id on otus (marker_wise_sequence_id)")

                c.execute("CREATE INDEX markers_marker on markers (marker)")

                logging.info("Creating SQL indexes on nucleotides ..")
                c.execute("CREATE INDEX nucleotides_marker_id on nucleotides (marker_id)")
                c.execute("CREATE INDEX nucleotides_sequence on nucleotides (sequence)")
                c.execute("CREATE INDEX nucleotides_marker_wise_id on nucleotides (marker_wise_id)")

                logging.info("Creating SQL indexes on taxonomy ..")
                c.execute("CREATE INDEX taxonomy_taxonomy on taxonomy (taxonomy)")
                db.commit()


                logging.info("Creating sorted protein sequences data ..")
                # Write a file of nucleotide id + protein sequence, and sort
                sorted_proteins_path = os.path.join(my_tempdir,'makedb_sort_output_protein')
                # Have to set LC_COLLATE=C because otherwise dashes in sequences can be ignored
                proc = subprocess.Popen(['bash','-c','LC_COLLATE=C sort --parallel={} --buffer-size=20% > {}'.format(num_threads, sorted_proteins_path)],
                    stdin=subprocess.PIPE,
                    stdout=None,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
                with open(numbered_sequences_file) as csvfile_in:
                    reader = csv.reader(csvfile_in, delimiter="\t")

                    for row in reader:
                        # id INTEGER PRIMARY KEY,"
                        # " marker_id int, sequence text, marker_wise_id int
                        nucleotide_id = row[0]
                        marker_id = row[1]
                        sequence = row[2]
                        print("\t".join([marker_id, nucleotides_to_protein(sequence), nucleotide_id]), file=proc.stdin)
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise Exception("Sort command returned non-zero exit status %i.\n"\
                        "STDERR was: %s" % (
                            proc.returncode, proc.stderr.read()))
                proc.stderr.close()

                logging.info("Creating protein sequence data for import ..")
                # Write a file with unique protein IDs, and the join table between nuc and prot
                proteins_file = os.path.join(my_tempdir,'makedb_numbered_proteins.tsv')
                nucleotide_proteins_file = os.path.join(my_tempdir,'nucleotides_proteins.tsv')
                with open(sorted_proteins_path) as csvfile_in:
                    reader = csv.reader(csvfile_in, delimiter="\t")
                    last_protein = None
                    last_protein_id = 0
                    last_marker_id = None
                    last_marker_wise_protein_id = None

                    with open(proteins_file,'w') as proteins_file_io:
                        with open(nucleotide_proteins_file,'w') as nucleotide_proteins_file_io:
                            for i, row in enumerate(reader):
                                marker_id = row[0]
                                current_protein_sequence = row[1]
                                nucleotide_id = row[2]
                                if current_protein_sequence != last_protein:
                                    last_protein = current_protein_sequence
                                    last_protein_id += 1
                                    if last_marker_id != marker_id:
                                        last_marker_wise_protein_id = 0 #start from 0 for annoy/scann
                                        last_marker_id = marker_id
                                    print("\t".join([str(last_protein_id), str(last_marker_wise_protein_id), current_protein_sequence]), file=proteins_file_io)
                                    last_marker_wise_protein_id += 1
                                print("\t".join([str(i), nucleotide_id, str(last_protein_id)]), file=nucleotide_proteins_file_io)
                logging.info("Running imports of proteins and protein_nucleotides ..")
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE proteins ("
                    "id INTEGER PRIMARY KEY AUTOINCREMENT,"
                    "marker_wise_id INTEGER, "
                    "protein_sequence text);\n"
                        '.separator "\\t"\n'
                        ".import {} proteins".format(proteins_file))
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE nucleotides_proteins ("
                    "id INTEGER PRIMARY KEY AUTOINCREMENT,"
                    "nucleotide_id int,"
                    "protein_id int);\n"
                        '.separator "\\t"\n'
                        ".import {} nucleotides_proteins".format(nucleotide_proteins_file))

                logging.info("Creating proteins indices ..")
                c.execute("CREATE INDEX proteins_sequence on proteins (protein_sequence)")
                c.execute("CREATE INDEX proteins_marker_wise_id on proteins (marker_wise_id)")
                db.commit()

                logging.info("Creating nucleotides_proteins indices ..")
                c.execute("CREATE INDEX nucleotides_proteins_protein_id on nucleotides_proteins (protein_id)")
                c.execute("CREATE INDEX nucleotides_proteins_nucleotide_id on nucleotides_proteins (nucleotide_id)")
                db.commit()
            

        # Create sequence indices
        sdb = SequenceDatabase.acquire(db_path)
        if 'scann-naive' in sequence_database_methods:
            sequence_database_methods.append(SCANN_NAIVE_INDEX_FORMAT)
        if SCANN_INDEX_FORMAT in sequence_database_methods or SCANN_NAIVE_INDEX_FORMAT in sequence_database_methods:
            sdb.create_scann_indexes(sequence_database_types, SCANN_NAIVE_INDEX_FORMAT in sequence_database_methods)

        if NMSLIB_INDEX_FORMAT in sequence_database_methods:
            if NUCLEOTIDE_DATABASE_TYPE in sequence_database_types:
                sdb.create_nmslib_nucleotide_indexes()
            if PROTEIN_DATABASE_TYPE in sequence_database_types:
                sdb.create_nmslib_protein_indexes()

        if ANNOY_INDEX_FORMAT in sequence_database_methods:
            if NUCLEOTIDE_DATABASE_TYPE in sequence_database_types:
                sdb.create_annoy_nucleotide_indexes(ntrees=num_annoy_nucleotide_trees)
            if PROTEIN_DATABASE_TYPE in sequence_database_types:
                sdb.create_annoy_protein_indexes(ntrees=num_annoy_protein_trees)

        if SMAFA_NAIVE_INDEX_FORMAT in sequence_database_methods:
            if NUCLEOTIDE_DATABASE_TYPE in sequence_database_types:
                sdb.create_smafa_naive_nucleotide_indexes()

        logging.info("Finished singlem DB creation")

    def create_smafa_naive_nucleotide_indexes(self):
        logging.info("Creating smafa-naive nucleotide sequence indices ..")
        nucleotide_db_dir = os.path.join(self.base_directory, 'nucleotide_indices_smafa_naive')
        os.makedirs(nucleotide_db_dir)
        
        for marker_row in self.sqlalchemy_connection.execute(select(Marker)):
            marker_name = marker_row.marker
            logging.info("Tabulating unique nucleotide sequences for {}..".format(marker_name))
            count = 0

            with tempfile.NamedTemporaryFile(prefix='singlem-smafa-create-', suffix='.fasta') as fasta_file:
                for row in self.sqlalchemy_connection.execute(select(
                    NucleotideSequence.sequence, NucleotideSequence.marker_wise_id) \
                    .where(NucleotideSequence.marker_id == marker_row.id)
                    .order_by(NucleotideSequence.marker_wise_id)):

                    fasta_file.write(str.encode(">{}\n{}\n".format(row.marker_wise_id, row.sequence)))
                    count += 1
                fasta_file.flush()

                extern.run('smafa makedb --database {} --input {}'.format(
                    os.path.join(nucleotide_db_dir, "%s.smafa_naive_index" % marker_name),
                    fasta_file.name))
            logging.info("Finished writing index containing {} sequences to disk".format(count))


    def create_nmslib_nucleotide_indexes(self):
        logging.info("Creating nmslib nucleotide sequence indices ..")
        nucleotide_db_dir = os.path.join(self.base_directory, 'nucleotide_indices_nmslib')
        os.makedirs(nucleotide_db_dir)
        
        for marker_row in self.sqlalchemy_connection.execute(select(Marker)):
            nucleotide_index = SequenceDatabase._nucleotide_nmslib_init()

            marker_name = marker_row.marker
            logging.info("Tabulating unique nucleotide sequences for {}..".format(marker_name))
            count = 0

            for row in self.sqlalchemy_connection.execute(select(
                NucleotideSequence.sequence, NucleotideSequence.marker_wise_id) \
                .where(NucleotideSequence.marker_id == marker_row.id)):

                nucleotide_index.addDataPoint(row.marker_wise_id, nucleotides_to_binary(row.sequence))
                count += 1

            # TODO: Tweak index creation parameters?
            logging.info("Creating binary nucleotide index from {} unique sequences ..".format(count))
            nucleotide_index.createIndex()

            logging.info("Writing index to disk ..")
            nucleotide_db_path = os.path.join(nucleotide_db_dir, "%s.nmslib_index" % marker_name)
            nucleotide_index.saveIndex(nucleotide_db_path, save_data=True)
            logging.info("Finished writing index to disk")

    def create_nmslib_protein_indexes(self):
        logging.info("Creating nmslib protein sequence indices ..")
        protein_db_dir = os.path.join(self.base_directory, 'protein_indices_nmslib')
        os.makedirs(protein_db_dir)
        for marker_row in self.sqlalchemy_connection.execute(select(Marker)):
            protein_index = SequenceDatabase._protein_nmslib_init()

            marker_name = marker_row.marker
            logging.info("Tabulating unique protein sequences for {}..".format(marker_name))
            count = 0

            for row in self.sqlalchemy_connection.execute(select(
                distinct(ProteinSequence.marker_wise_id), ProteinSequence.protein_sequence) \
                    .where(ProteinSequence.id == NucleotidesProteins.protein_id) \
                    .where(NucleotidesProteins.nucleotide_id == NucleotideSequence.id) \
                    .where(NucleotideSequence.marker_id == marker_row.id)):
                protein_index.addDataPoint(row.marker_wise_id, protein_to_binary(row.protein_sequence))
                count += 1

            # TODO: Tweak index creation parameters?
            logging.info("Creating binary protein index from {} unique sequences ..".format(count))
            protein_index.createIndex()

            logging.info("Writing index to disk ..")
            protein_db_path = os.path.join(protein_db_dir, "%s.nmslib_index" % marker_name)
            protein_index.saveIndex(protein_db_path, save_data=True)
            logging.info("Finished writing index to disk")

    def create_annoy_nucleotide_indexes(self, ntrees):
        logging.info("Creating annoy nucleotide sequence indices ..")
        nucleotide_db_dir = os.path.join(self.base_directory, 'nucleotide_indices_annoy')
        os.makedirs(nucleotide_db_dir)

        for marker_row in self.sqlalchemy_connection.execute(select(Marker)):
            annoy_index = self._nucleotide_annoy_init()

            marker_name = marker_row.marker
            logging.info("Tabulating unique nucleotide sequences for {}..".format(marker_name))
            count = 0

            for row in self.sqlalchemy_connection.execute(select(
                NucleotideSequence.sequence, NucleotideSequence.marker_wise_id) \
                .where(NucleotideSequence.marker_id == marker_row.id)):

                annoy_index.add_item(row.marker_wise_id, nucleotides_to_binary_array(row.sequence))
                count += 1

            # TODO: Tweak index creation parameters?
            logging.info("Creating binary nucleotide index from {} unique sequences and ntrees={}..".format(count, ntrees))
            annoy_index.build(ntrees)

            logging.info("Writing index to disk ..")
            annoy_index.save(os.path.join(nucleotide_db_dir, "%s.annoy_index" % marker_name))
            logging.info("Finished writing index to disk")
            # Delete immediately to save RAM (was using 200G+ before getting killed on big DB)
            del annoy_index

    def create_annoy_protein_indexes(self, ntrees):
        logging.info("Creating annoy protein sequence indices ..")
        protein_db_dir = os.path.join(self.base_directory, 'protein_indices_annoy')
        os.makedirs(protein_db_dir)

        for marker_row in self.sqlalchemy_connection.execute(select(Marker)):
            annoy_index = self._protein_annoy_init()

            marker_name = marker_row.marker
            logging.info("Tabulating unique protein sequences for {}..".format(marker_name))
            count = 0

            for row in self.sqlalchemy_connection.execute(select(
                distinct(ProteinSequence.marker_wise_id), ProteinSequence.protein_sequence) \
                    .where(ProteinSequence.id == NucleotidesProteins.protein_id) \
                    .where(NucleotidesProteins.nucleotide_id == NucleotideSequence.id) \
                    .where(NucleotideSequence.marker_id == marker_row.id)):

                annoy_index.add_item(row.marker_wise_id, protein_to_binary_array(row.protein_sequence))
                count += 1

            # TODO: Tweak index creation parameters?
            logging.info("Creating binary protein index from {} unique sequences and ntrees={}..".format(count, ntrees))
            annoy_index.build(ntrees)

            logging.info("Writing index to disk ..")
            annoy_index.save(os.path.join(protein_db_dir, "%s.annoy_index" % marker_name))
            logging.info("Finished writing index to disk")
            # Delete immediately to save RAM (was using 200G+ before getting killed on big DB)
            del annoy_index

    def create_scann_indexes(self, sequence_database_types, generate_brute_force_index):
        logging.info("Creating scann sequence indices ..")
        import tensorflow as tf
        import scann # only load when needed to speed start-up
        if NUCLEOTIDE_DATABASE_TYPE in sequence_database_types:
            nucleotide_db_dir_ah = os.path.join(self.base_directory, 'nucleotide_indices_scann')
            nucleotide_db_dir_brute_force = os.path.join(self.base_directory, 'nucleotide_indices_scann_brute_force')
            os.makedirs(nucleotide_db_dir_ah)
            if generate_brute_force_index:
                os.makedirs(nucleotide_db_dir_brute_force)
        if PROTEIN_DATABASE_TYPE in sequence_database_types:
            protein_db_dir_ah = os.path.join(self.base_directory, 'protein_indices_scann')
            protein_db_dir_brute_force = os.path.join(self.base_directory, 'protein_indices_scann_brute_force')
            os.makedirs(protein_db_dir_ah)
            if generate_brute_force_index:
                os.makedirs(protein_db_dir_brute_force)

        def generate_indices_from_array(a, db_dir_ah, db_dir_brute_force, generate_brute_force_index):
            normalized_dataset = a / np.linalg.norm(a, axis=1)[:, np.newaxis]
            logging.info("Found {} sequences for {}".format(a.shape[0], marker_name))
            del a # not sure if this matters much

            logging.info("Creating SCANN AH index ..")
            searcher = scann.scann_ops_pybind.builder(normalized_dataset, 10, "dot_product").tree(
                num_leaves=round(np.sqrt(normalized_dataset.shape[0])), num_leaves_to_search=100, training_sample_size=250000).score_ah(
                2, anisotropic_quantization_threshold=0.2).reorder(100).build()
            directory = os.path.join(db_dir_ah, marker_name)
            os.mkdir(directory)
            searcher.serialize(directory)
            del searcher

            if generate_brute_force_index:
                logging.info("Creating SCANN brute force index ..")
                # use scann.scann_ops.build() to instead create a
                # TensorFlow-compatible searcher could not work out how to
                # deserialise a brute force without any doco, so just copying method
                # from the tests i.e.
                # https://github.com/google-research/google-research/blob/34444253e9f57cd03364bc4e50057a5abe9bcf17/scann/scann/scann_ops/py/scann_ops_test.py#L93
                searcher_naive = scann.scann_ops.builder(normalized_dataset, 10, "dot_product").tree(
                    num_leaves=round(np.sqrt(normalized_dataset.shape[0])), num_leaves_to_search=100).score_brute_force(True).build()
                directory = os.path.join(db_dir_brute_force, marker_name)
                module = searcher_naive.serialize_to_module()
                tf.saved_model.save(
                    module,
                    directory,
                    options=tf.saved_model.SaveOptions(namespace_whitelist=["Scann"]))
                del searcher_naive

        for marker_row in self.sqlalchemy_connection.execute(select(Marker)):
            marker_name = marker_row.marker
            marker_id = marker_row.id
            
            if NUCLEOTIDE_DATABASE_TYPE in sequence_database_types:
                logging.info("Tabulating unique nucleotide sequences for {}..".format(marker_name))
                a = np.concatenate([np.array([nucleotides_to_binary_array(entry.sequence)]) for entry in \
                    self.sqlalchemy_connection.execute(select(
                        NucleotideSequence.sequence) \
                        .where(NucleotideSequence.marker_id == marker_id) \
                        .order_by(NucleotideSequence.marker_wise_id))
                ])
                if a.shape[0] < 16:
                    logging.warning("Adding dummy nucleotide sequences to SCANN AH/NAIVE DB creation since the number of real datapoints is too small")
                    a = np.concatenate([a, np.ones((16-a.shape[0], a.shape[1]))])
                generate_indices_from_array(a, nucleotide_db_dir_ah, nucleotide_db_dir_brute_force, generate_brute_force_index)
                logging.info("Finished writing nucleotide indices to disk")
            
            if PROTEIN_DATABASE_TYPE in sequence_database_types:
                logging.info("Tabulating unique protein sequences for {}..".format(marker_name))
                a = np.concatenate([np.array([protein_to_binary_array(entry.protein_sequence)]) for entry in \
                    self.sqlalchemy_connection.execute(select(
                        ProteinSequence.protein_sequence) \
                            .order_by(ProteinSequence.marker_wise_id) \
                            .where(ProteinSequence.id == NucleotidesProteins.protein_id) \
                            .where(NucleotidesProteins.nucleotide_id == NucleotideSequence.id) \
                            .where(NucleotideSequence.marker_id == marker_id)
                            .distinct())
                ])
                if a.shape[0] < 16:
                    logging.warn("Adding dummy protein sequences to SCANN AH/NAIVE DB creation since the number of real datapoints is too small")
                    a = np.concatenate([a, np.ones((16-a.shape[0], a.shape[1]))])
                generate_indices_from_array(a, protein_db_dir_ah, protein_db_dir_brute_force, generate_brute_force_index)
                del a
                logging.info("Finished writing protein indices to disk")

    
    @staticmethod
    def dump(db_path):
        """Dump the DB contents to STDOUT, requiring a version 5+ database"""
        db = SequenceDatabase.acquire(db_path)
        table = SequenceDatabaseOtuTable(db)
        
        print("\t".join(OtuTable.DEFAULT_OUTPUT_FIELDS))
        # DEFAULT_OUTPUT_FIELDS = str.split('gene sample sequence num_hits coverage taxonomy')

        for otu in table:
            print("\t".join([
                otu.marker,
                otu.sample_name,
                otu.sequence,
                str(otu.count),
                '%.2f' % otu.coverage,
                otu.taxonomy
            ]))


class SequenceDatabaseOtuTable:

    def __init__(self, db):
        self.db = db

    def __iter__(self):
        engine = self.db.engine

        with engine.connect() as conn:
            # First cache the taxonomy
            taxonomy_entries = Taxonomy.generate_python_index(conn)

            # And markers
            marker_entries = Marker.generate_python_index(conn)
        
            batch_size = 10000
            builder = select(
                Otu.marker_id, Otu.sample_name, Otu.sequence, Otu.num_hits, Otu.coverage, Otu.taxonomy_id
                ).execution_options(yield_per=batch_size)
        
            for batch in conn.execute(builder).partitions(batch_size):
                for row in batch:
                    entry = OtuTableEntry()
                    # marker = None
                    # sample_name = None
                    # sequence = None
                    # count = None
                    # taxonomy = None
                    # coverage = None
                    # data = None
                    # fields = None
                    entry.marker = marker_entries[row.marker_id]
                    entry.sample_name = row.sample_name
                    entry.sequence = row.sequence
                    entry.count = row.num_hits
                    entry.coverage = row.coverage
                    entry.taxonomy = taxonomy_entries[row.taxonomy_id]
                    yield entry


def _base_to_binary(x):
    if x == 'A':
        return '1 0 0 0 0'
    elif x == 'T':
        return '0 1 0 0 0'
    elif x == 'C':
        return '0 0 1 0 0'
    elif x == 'G':
        return '0 0 0 1 0'
    else:
        return '0 0 0 0 1'
    
def nucleotides_to_binary(seq):
    return ' '.join([_base_to_binary(b) for b in seq])
    
def _base_to_binary_array(x):
    if x == 'A':
        return [1,0,0,0,0]
    elif x == 'T':
        return [0,1,0,0,0]
    elif x == 'C':
        return [0,0,1,0,0]
    elif x == 'G':
        return [0,0,0,1,0]
    else:
        return [0,0,0,0,1]
    
def nucleotides_to_binary_array(seq):
    return list(itertools.chain(*[_base_to_binary_array(b) for b in seq]))

AA_ORDER = ['W',
        'H',
        'Q',
        'M',
        'K',
        'G',
        'A',
        'I',
        'F',
        'S',
        'L',
        'Y',
        'T',
        'D',
        'E',
        'R',
        'P',
        'V',
        'N',
        'C',

        '-',
        'X']

def _aa_to_binary(x):
    return ' '.join(['1' if aa == x else '0' for aa in AA_ORDER])

def protein_to_binary(seq):
    return ' '.join([_aa_to_binary(b) for b in seq])
    
def _aa_to_binary_array(x):
    return list([1 if aa == x else 0 for aa in AA_ORDER])

def protein_to_binary_array(seq):
    return list(itertools.chain(*[_aa_to_binary_array(b) for b in seq]))

def nucleotides_to_protein(seq):
    aas = []
    codon_table=Bio.Data.CodonTable.standard_dna_table.forward_table
    for i in range(0, len(seq), 3):
        codon = seq[i:(i+3)]
        if codon == '---':
            aas.append('-')
        elif codon in codon_table:
            aas.append(codon_table[codon])
        else:
            aas.append('X')
    return ''.join(aas)
