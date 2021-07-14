import os
import tempfile
import logging
import subprocess
import sqlite3
import glob
import json
import itertools
import sys
import csv
import extern

# Import masonite, pretty clunky
sys.path = [os.path.dirname(os.path.realpath(__file__))] + sys.path
import config.database
from masoniteorm.connections import ConnectionResolver
from masoniteorm.query import QueryBuilder

import nmslib
import numba
from numba.core import types
import Bio.Data.CodonTable
from annoy import AnnoyIndex

from .otu_table import OtuTableEntry, OtuTable

DEFAULT_NUM_THREADS = 1


class SequenceDatabase:
    version = 4
    SQLITE_DB_NAME = 'otus.sqlite3'
    _marker_to_nmslib_nucleotide_index_file = {}
    _marker_to_nmslib_protein_index_file = {}
    _marker_to_annoy_nucleotide_index_file = {}
    _marker_to_annoy_protein_index_file = {}

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    VERSION_KEY = 'singlem_database_version'

    _REQUIRED_KEYS = {4: [VERSION_KEY]}

    def add_sequence_db(self, marker_name, db_path, index_format, sequence_type):
        if index_format == 'nmslib':
            if sequence_type == 'nucleotide':
                self._marker_to_nmslib_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == 'protein':
                self._marker_to_nmslib_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        elif index_format == 'annoy':
            if sequence_type == 'nucleotide':
                self._marker_to_annoy_nucleotide_index_file[marker_name] = db_path
            elif sequence_type == 'protein':
                self._marker_to_annoy_protein_index_file[marker_name] = db_path
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        else:
            raise Exception("Unknown index type {}".format(index_format))

    def get_sequence_index(self, marker_name, index_format, sequence_type):
        index = None
        if index_format == 'nmslib':
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
        elif index_format == 'annoy':
            if sequence_type == 'nucleotide':
                raise NotImplementedException()
            elif sequence_type == 'protein':
                raise NotImplementedException()
            else:
                raise Exception('Invalid sequence type: %s' % sequence_type)
        else:
            raise Exception("Unknown index type {}".format(index_format))
        
        logging.warn("No %s / %s sequence index DB found for %s" % (index_format, sequence_type, marker_name))
        return None

    @staticmethod
    def _nucleotide_nmslib_init():
        return nmslib.init(space='bit_hamming', data_type=nmslib.DataType.OBJECT_AS_STRING, dtype=nmslib.DistType.INT, method='hnsw')

    def query_builder(self, check=False):
        return SequenceDatabase._query_builder(self.sqlite_file, check=check)

    @staticmethod
    def _query_builder(sqlite_db_path, check=False):
        '''setup config and return a query builder. If check==True, then make sure there the connection is alive, raising an exception otherwise'''
        if check:
            logging.debug("Connecting to DB {}".format(sqlite_db_path))        
            if not os.path.exists(sqlite_db_path):
                raise Exception("SQLite3 database does not appear to exist in the SingleM database - perhaps it is the wrong version?")

        DB = ConnectionResolver().set_connection_details({
            'default': 'singlem',
            'singlem': {
                'driver': 'sqlite',
                'database': sqlite_db_path,
                'prefix': '',
                'log_queries': False
            }
        })
        config.database.DB = DB

        qb = QueryBuilder(connection='singlem', connection_details=DB.get_connection_details())

        if check:
            try:
                len(qb.table('otus').limit(1).get())
                qb = QueryBuilder(connection='singlem', connection_details=DB.get_connection_details())
            except Exception as e:
                logging.error("Failure to extract any data from the otus table of the SQ Lite DB indicates this SingleM DB is either too old or is corrupt.")
                raise(e)

        return qb

    @staticmethod
    def acquire(path):
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
        if found_version == 4:
            for key in SequenceDatabase._REQUIRED_KEYS[found_version]:
                if key not in db._contents_hash:
                    raise Exception(
                        "Unexpectedly did not find required key {} in SingleM database contents file: {}".format(
                            key, path))
        else:
            raise Exception("Unexpected SingleM DB version found: {}".format(found_version))

        db.sqlite_file = os.path.join(path, SequenceDatabase.SQLITE_DB_NAME)

        nmslib_nucleotide_index_files = glob.glob("%s/nucleotide_indices_nmslib/*.nmslib_index" % path)
        logging.debug("Found nmslib_nucleotide_index_files: %s" % ", ".join(nmslib_nucleotide_index_files))
        for g in nmslib_nucleotide_index_files:
            marker = os.path.basename(g).replace('.nmslib_index','')
            db.add_sequence_db(marker, g, 'nmslib','nucleotide')

        nmslib_protein_index_files = glob.glob("%s/protein_indices_nmslib/*.nmslib_index" % path)
        logging.debug("Found nmslib_protein_index_files: %s" % ", ".join(nmslib_protein_index_files))
        for g in nmslib_protein_index_files:
            marker = os.path.basename(g).replace('.nmslib_index','')
            db.add_sequence_db(marker, g, 'nmslib','protein')

        annoy_nucleotide_index_files = glob.glob("%s/nucleotide_indices_annoy/*.annoy_index" % path)
        logging.debug("Found annoy_nucleotide_index_files: %s" % ", ".join(annoy_nucleotide_index_files))
        for g in annoy_nucleotide_index_files:
            marker = os.path.basename(g).replace('.annoy_index','')
            db.add_sequence_db(marker, g, 'annoy','nucleotide')

        annoy_protein_index_files = glob.glob("%s/protein_indices_annoy/*.annoy_index" % path)
        logging.debug("Found annoy_protein_index_files: %s" % ", ".join(nmslib_protein_index_files))
        for g in annoy_protein_index_files:
            marker = os.path.basename(g).replace('.annoy_index','')
            db.add_sequence_db(marker, g, 'annoy','protein')
        
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
        tmpdir=None):

        if num_threads is None:
            num_threads = DEFAULT_NUM_THREADS
        if pregenerated_sqlite3_db:
            logging.info("Re-using previous SQLite database {}".format(pregenerated_sqlite3_db))
            sqlite_db_path = pregenerated_sqlite3_db

            marker_list = set()
            for row in SequenceDatabase._query_builder(sqlite_db_path).table('otus').select_raw("distinct(marker) as marker").get():
                marker_list.add(row['marker'])
            logging.info("Found {} markers e.g. {}".format(len(marker_list), list(marker_list)[0]))

        else:
            # ensure db does not already exist
            if os.path.exists(db_path):
                raise Exception("Cowardly refusing to overwrite already-existing database path '%s'" % db_path)
            logging.info("Creating SingleM database at {}".format(db_path))
            os.makedirs(db_path)

            # Create contents file
            contents_file_path = os.path.join(db_path, SequenceDatabase._CONTENTS_FILE_NAME)
            with open(contents_file_path, 'w') as f:
                json.dump({
                    SequenceDatabase.VERSION_KEY: 4,
                }, f)

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
                # create tempdir
                sorted_path = os.path.join(my_tempdir,'makedb_sort_output')
                proc = subprocess.Popen(['bash','-c','sort --parallel={} --buffer-size=20% > {}'.format(num_threads, sorted_path)],
                    stdin=subprocess.PIPE,
                    stdout=None,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
                for entry in otu_table_collection:
                    total_otu_count += 1
                    print("\t".join((entry.marker, entry.sequence, entry.sample_name, str(entry.count),
                                str(entry.coverage), entry.taxonomy)), file=proc.stdin)
                logging.info("Sorting {} OTU observations ..".format(total_otu_count))
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise Exception("Sort command returned non-zero exit status %i.\n"\
                        "STDERR was: %s" % (
                            proc.returncode, proc.stderr.read()))

                logging.info("Creating numbered OTU table tsv")
                marker_index = 0
                sequence_index = 0
                numbered_table_file = os.path.join(my_tempdir,'makedb_numbered_output.tsv')
                numbered_marker_and_sequence_file = os.path.join(my_tempdir,'number_and_sequence_file')
                with open(sorted_path) as csvfile_in:
                    reader = csv.reader(csvfile_in, delimiter="\t")
                    last_marker = None
                    last_sequence = None
                    with open(numbered_table_file,'w') as otus_output_table_io:
                        with open(numbered_marker_and_sequence_file,'w') as marker_and_sequence_foutput_table_io:
                            for i, row in enumerate(reader):
                                if last_marker != row[0]:
                                    last_marker = row[0]
                                    marker_index += 1
                                    last_sequence = row[1]
                                    sequence_index += 1
                                elif last_sequence != row[1]:
                                    last_sequence = row[1]
                                    sequence_index += 1
                                print("\t".join([str(i+1)]+row[2:]+[str(marker_index),str(sequence_index)]), file=otus_output_table_io)
                                print("\t".join([str(marker_index),str(sequence_index),row[0],row[1]]), file=marker_and_sequence_foutput_table_io)

                logging.info("Importing OTU table into SQLite ..")
                sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE otus (id INTEGER PRIMARY KEY AUTOINCREMENT,"
                        " sample_name text, num_hits int, coverage float, taxonomy text, marker_id int, sequence_id int);\n"
                        '.separator "\\t"\n'
                        ".import {} otus".format(numbered_table_file))

                logging.info("Creating markers and nucleotide table TSV files ..")
                numbered_markers_file = os.path.join(my_tempdir,'makedb_numbered_markers_file')
                numbered_sequences_file = os.path.join(my_tempdir,'makedb_numbered_sequences_file')
                with open(numbered_marker_and_sequence_file) as csvfile_in:
                    reader = csv.reader(csvfile_in, delimiter="\t")
                    last_marker_id = None
                    last_sequence_id = None
                    last_marker_wise_sequence_id = None
                    with open(numbered_markers_file,'w') as markers_output_table_io:
                        with open(numbered_sequences_file,'w') as nucleotides_output_table_io:
                            for row in reader:
                                if last_marker_id != row[0]:
                                    last_marker_wise_sequence_id = 1
                                    print("\t".join([row[0], row[2]]), file=markers_output_table_io)
                                    print("\t".join([row[1], row[0], row[3], str(last_marker_wise_sequence_id)]), file=nucleotides_output_table_io)
                                    last_marker_id = row[0]
                                    last_sequence_id = row[1]
                                elif last_sequence_id != row[1]:
                                    print("\t".join([row[1], row[0], row[3], str(last_marker_wise_sequence_id)]), file=nucleotides_output_table_io)
                                    last_sequence_id = row[1]
                                last_marker_wise_sequence_id += 1

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
                c.execute("CREATE INDEX otu_taxonomy on otus (taxonomy)")
                c.execute("CREATE INDEX otu_marker on otus (marker_id)")
                c.execute("CREATE INDEX otu_sequence on otus (sequence_id)")

                c.execute("CREATE INDEX markers_marker on markers (marker)")

                logging.info("Creating SQL indexes on nucleotides ..")
                c.execute("CREATE INDEX nucleotides_marker_id on nucleotides (marker_id)")
                c.execute("CREATE INDEX nucleotides_sequence on nucleotides (sequence)")
                db.commit()


                logging.info("Creating sorted protein sequences data ..")
                # Write a file of nucleotide id + protein sequence, and sort
                sorted_proteins_path = os.path.join(my_tempdir,'makedb_sort_output_protein')
                proc = subprocess.Popen(['bash','-c','sort --parallel={} --buffer-size=20% > {}'.format(num_threads, sorted_proteins_path)],
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
                        sequence = row[2]
                        print("\t".join([nucleotides_to_protein(row[2]), nucleotide_id]), file=proc.stdin)
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise Exception("Sort command returned non-zero exit status %i.\n"\
                        "STDERR was: %s" % (
                            proc.returncode, proc.stderr.read()))

                logging.info("Creating protein sequence data for import ..")
                # Write a file with unique protein IDs, and the join table between nuc and prot
                proteins_file = os.path.join(my_tempdir,'makedb_numbered_proteins.tsv')
                nucleotide_proteins_file = os.path.join(my_tempdir,'nucleotides_proteins.tsv')
                with open(sorted_proteins_path) as csvfile_in:
                    reader = csv.reader(csvfile_in, delimiter="\t")
                    last_protein = None
                    last_protein_id = 0

                    with open(proteins_file,'w') as proteins_file_io:
                        with open(nucleotide_proteins_file,'w') as nucleotide_proteins_file_io:
                            for i, row in enumerate(reader):
                                current_protein_sequence = row[0]
                                nucleotide_id = row[1]
                                if current_protein_sequence != last_protein:
                                    last_protein = current_protein_sequence
                                    last_protein_id += 1
                                    print("\t".join([str(last_protein_id), current_protein_sequence]), file=proteins_file_io)
                                print("\t".join([str(i), str(last_protein_id), nucleotide_id]), file=nucleotide_proteins_file_io)
                logging.info("Running imports of proteins and protein_nucleotides ..")
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE proteins ("
                    "id INTEGER PRIMARY KEY AUTOINCREMENT,"
                    "protein_sequence int);\n"
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
                db.commit()

                logging.info("Creating nucleotides_proteins indices ..")
                c.execute("CREATE INDEX nucleotides_proteins_protein_id on nucleotides_proteins (protein_id)")
                c.execute("CREATE INDEX nucleotides_proteins_nucleotide_id on nucleotides_proteins (nucleotide_id)")
                db.commit()
            

        # Create nucleotide index files
        SequenceDatabase.acquire(db_path).create_nmslib_nucleotide_indexes()
        SequenceDatabase.acquire(db_path).create_nmslib_protein_indexes()

        SequenceDatabase.acquire(db_path).create_annoy_nucleotide_indexes(ntrees=20)
        SequenceDatabase.acquire(db_path).create_annoy_protein_indexes(ntrees=20)

        logging.info("Finished singlem DB creation")

    def create_nmslib_nucleotide_indexes(self):
        logging.info("Creating nmslib nucleotide sequence indices ..")
        nucleotide_db_dir = os.path.join(self.base_directory, 'nucleotide_indices_nmslib')
        os.makedirs(nucleotide_db_dir)
        
        for marker_row in self.query_builder().table('markers').get():
            nucleotide_index = SequenceDatabase._nucleotide_nmslib_init()

            marker_name = marker_row['marker']
            logging.info("Tabulating unique nucleotide sequences for {}..".format(marker_name))
            count = 0

            for row in self.query_builder().table('nucleotides').select('sequence').select('id').where('marker_id', marker_row['id']).get():
                # logging.debug("Adding sequence with ID {}: {}".format(row['id'], row['sequence']))
                nucleotide_index.addDataPoint(row['id'], nucleotides_to_binary(row['sequence']))
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
        for marker_row in self.query_builder().table('markers').get():
            protein_index = SequenceDatabase._nucleotide_nmslib_init()

            marker_name = marker_row['marker']
            logging.info("Tabulating unique protein sequences for {}..".format(marker_name))
            count = 0

            for row in self.query_builder().table('proteins'). \
                select_raw('proteins.id as id, protein_sequence').join('nucleotides_proteins','proteins.id','=','nucleotides_proteins.protein_id'). \
                    join('nucleotides','nucleotides_proteins.nucleotide_id','=','nucleotides.id'). \
                        where('nucleotides.marker_id', marker_row['id']).get():
                protein_index.addDataPoint(row['id'], protein_to_binary(row['protein_sequence']))
                count += 1

            # TODO: Tweak index creation parameters?
            logging.info("Creating binary protein index from {} unique sequences ..".format(count))
            protein_index.createIndex()

            logging.info("Writing index to disk ..")
            protein_db_path = os.path.join(protein_db_dir, "%s.nmslib_index" % marker_name)
            protein_index.saveIndex(protein_db_path, save_data=True)
            logging.info("Finished writing index to disk")

    def create_annoy_nucleotide_indexes(self, ntrees=None):
        example_seq = self.query_builder().table('nucleotides').limit(1).first()['sequence']
        ndim = len(example_seq)*5
        logging.debug("Creating {} dimensional annoy indices ..".format(ndim))

        logging.info("Creating annoy nucleotide sequence indices ..")
        nucleotide_db_dir = os.path.join(self.base_directory, 'nucleotide_indices_annoy')
        os.makedirs(nucleotide_db_dir)

        for marker_row in self.query_builder().table('markers').get():
            annoy_index = AnnoyIndex(ndim, 'hamming')

            marker_name = marker_row['marker']
            logging.info("Tabulating unique nucleotide sequences for {}..".format(marker_name))
            count = 0

            for row in self.query_builder().table('nucleotides').select('sequence').select('marker_wise_id').where('marker_id', marker_row['id']).get():
                logging.debug("Adding sequence with ID {}: {}".format(row['marker_wise_id'], row['sequence']))
                annoy_index.add_item(row['marker_wise_id'], nucleotides_to_binary_array(row['sequence']))
                count += 1

            # TODO: Tweak index creation parameters?
            if ntrees is None:
                ntrees = int(count / 10) #complete guess atm
                if ntrees < 1:
                    ntrees = 1
            logging.info("Creating binary nucleotide index from {} unique sequences and ntrees={}..".format(count, ntrees))
            annoy_index.build(ntrees)

            logging.info("Writing index to disk ..")
            annoy_index.save(os.path.join(nucleotide_db_dir, "%s.annoy_index" % marker_name))
            logging.info("Finished writing index to disk")

    def create_annoy_protein_indexes(self, ntrees=None):
        example_seq = self.query_builder().table('proteins').limit(1).first()['protein_sequence']
        ndim = len(example_seq)*len(AA_ORDER)
        logging.debug("Creating {} dimensional annoy indices ..".format(ndim))

        logging.info("Creating annoy protein sequence indices ..")
        protein_db_dir = os.path.join(self.base_directory, 'protein_indices_annoy')
        os.makedirs(protein_db_dir)

        for marker_row in self.query_builder().table('markers').get():
            annoy_index = AnnoyIndex(ndim, 'hamming')

            marker_name = marker_row['marker']
            logging.info("Tabulating unique protein sequences for {}..".format(marker_name))
            count = 0

            for row in self.query_builder().table('proteins'). \
                join('nucleotides_proteins','proteins.id','=','nucleotides_proteins.protein_id'). \
                join('nucleotides','nucleotides_proteins.nucleotide_id','=','nucleotides.id'). \
                select('protein_sequence').select_raw('nucleotides.marker_wise_id as marker_wise_id'). \
                where('nucleotides.marker_id', marker_row['id']).get():
                logging.debug("Adding sequence with ID {}: {}".format(row['marker_wise_id'], row['protein_sequence']))
                annoy_index.add_item(row['marker_wise_id'], protein_to_binary_array(row['protein_sequence']))
                count += 1

            # TODO: Tweak index creation parameters?
            if ntrees is None:
                ntrees = int(count / 10) #complete guess atm
                if ntrees < 1:
                    ntrees = 1
            logging.info("Creating binary protein index from {} unique sequences and ntrees={}..".format(count, ntrees))
            annoy_index.build(ntrees)

            logging.info("Writing index to disk ..")
            annoy_index.save(os.path.join(protein_db_dir, "%s.annoy_index" % marker_name))
            logging.info("Finished writing index to disk")
    
    @staticmethod
    def dump(db_path):
        """Dump the DB contents to STDOUT, requiring only that the DB is a version that
        has an otus table in sqlite3 form (i.e. version 2 at least).

        """
        sqlite_db_path = os.path.join(db_path, SequenceDatabase.SQLITE_DB_NAME)
        
        print("\t".join(OtuTable.DEFAULT_OUTPUT_FIELDS))
        for chunk in SequenceDatabase._query_builder(sqlite_db_path).table('otus').chunk(1000):
            for entry in chunk:
                otu = OtuTableEntry()
                otu.marker = entry['marker']
                otu.sample_name = entry['sample_name']
                otu.sequence = entry['sequence']
                otu.count = entry['num_hits']
                otu.coverage = entry['coverage']
                otu.taxonomy = entry['taxonomy']
                print(str(otu))

@numba.njit()
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
    
@numba.njit()
def nucleotides_to_binary(seq):
    return ' '.join([_base_to_binary(b) for b in seq])
    
@numba.njit()
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
    
# @numba.njit()
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

# @numba.njit() # would like to do this, but better to move to lists not dict for codon table
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
