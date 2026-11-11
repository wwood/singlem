import logging
import re
import tempfile

from bird_tool_utils import iterable_chunks

from .singlem_package import SingleMPackage
from .biolib_lite.common import remove_extension

_GENOME_ACCESSION_REGEX = re.compile(r'(GC[AF]_\d+\.\d+)')
# supplement names its reference sequences '<spkg prefix>~<genome>‡<gene>'. The
# prefix never contains '~', but the genome name may.
_SUPPLEMENTED_GENOME_REGEX = re.compile(r'^[^~]*~(.+?)‡')


def genome_key_from_read_name(read_name):
    '''The key identifying the genome a metapackage reference sequence came from,
    or None. It is the sequence's versioned GCA_/GCF_ accession if it has one
    (e.g. 'rplB~GB_GCA_000013625.1' -> 'GCA_000013625.1'). Otherwise, for a
    sequence added by 'singlem supplement', it is the genome name supplement
    recorded: the FASTA basename with its last extension removed
    ('rplB~AMXMAG_0092.fa‡gene' -> 'AMXMAG_0092.fa'). Must agree with
    genome_key_from_genome_file.'''
    match = _GENOME_ACCESSION_REGEX.search(read_name)
    if match:
        return match.group(1)
    match = _SUPPLEMENTED_GENOME_REGEX.search(read_name)
    if match:
        return match.group(1)
    return None


def genome_key_from_genome_file(genome_file):
    '''The key identifying a genome FASTA file, e.g. a weebill profile's
    Genome_file. It is the versioned GCA_/GCF_ accession in the path if there is
    one. Otherwise it is the basename with its last extension removed, which is
    how 'singlem supplement' names genomes ('/x/AMXMAG_0092.fa.gz' ->
    'AMXMAG_0092.fa'). Must agree with genome_key_from_read_name.'''
    match = _GENOME_ACCESSION_REGEX.search(genome_file)
    if match:
        return match.group(1)
    return remove_extension(genome_file)

READ_NAME_TAXONOMY_TABLE = 'read_name_taxonomy'
TAXONOMY_MARKER_COUNT_TABLE = 'taxonomy_marker_count'


class MetapackageReadNameStore:
    '''Maps marker reference sequence names (and, via the embedded genome
    accession, genomes) to GTDB taxonomy. Stored as a DuckDB database in v7+
    metapackages and a SQLite database in earlier ones; the backend is detected
    from the database file's extension so old metapackages remain readable.'''

    def __init__(self, db_path, backend):
        self._db_path = db_path
        self._backend = backend

    @staticmethod
    def acquire(db_path):
        backend = 'duckdb' if str(db_path).endswith('.duckdb') else 'sqlite'
        return MetapackageReadNameStore(db_path, backend)

    def _connect(self):
        if self._backend == 'duckdb':
            import duckdb
            return duckdb.connect(self._db_path, read_only=True)
        import sqlite3
        return sqlite3.connect(self._db_path)

    @staticmethod
    def generate(singlem_package_paths, db_path, taxonomy_marker_counts=None):
        '''Build the read-name taxonomy store as a DuckDB database at db_path.'''
        import duckdb

        num_packages = 0
        num_read_names = 0
        with tempfile.NamedTemporaryFile(prefix='singlem_metapackage_read_name_store', suffix='.tsv') as temp_file:
            for singlem_package_path in singlem_package_paths:
                singlem_package = SingleMPackage.acquire(singlem_package_path)
                num_packages += 1
                for name, taxonomy in singlem_package.taxonomy_hash().items():
                    num_read_names += 1
                    temp_file.write("{}\t{}\t{}\n".format(
                        num_read_names, name, ';'.join(taxonomy)).encode('utf-8'))
            temp_file.flush()

            con = duckdb.connect(db_path)
            try:
                con.execute(
                    "CREATE TABLE {} (id BIGINT, read_name VARCHAR, taxonomy VARCHAR)".format(READ_NAME_TAXONOMY_TABLE))
                con.execute(
                    "COPY {} FROM '{}' (DELIMITER '\t', HEADER false, QUOTE '')".format(
                        READ_NAME_TAXONOMY_TABLE, temp_file.name))

                if taxonomy_marker_counts is not None:
                    logging.debug("Creating taxonomy_marker_count table")
                    with tempfile.NamedTemporaryFile(prefix='singlem_metapackage_read_name_store', suffix='.tsv') as temp_file2:
                        num_species = 0
                        for taxonomy, marker_count in taxonomy_marker_counts.items():
                            num_species += 1
                            temp_file2.write("{}\t{}\t{}\n".format(
                                num_species, taxonomy, marker_count).encode('utf-8'))
                        temp_file2.flush()
                        con.execute(
                            "CREATE TABLE {} (id BIGINT, taxonomy VARCHAR, marker_count BIGINT)".format(TAXONOMY_MARKER_COUNT_TABLE))
                        con.execute(
                            "COPY {} FROM '{}' (DELIMITER '\t', HEADER false, QUOTE '')".format(
                                TAXONOMY_MARKER_COUNT_TABLE, temp_file2.name))
            finally:
                con.close()
        logging.info("Imported {} packages and {} read names.".format(num_packages, num_read_names))

    def get_taxonomy_of_reads(self, read_names):
        '''Return dict of read name to taxonomy string'''
        to_return = {}
        conn = self._connect()
        try:
            for read_name_set in iterable_chunks(read_names, 900):
                names = [r for r in read_name_set if r is not None]
                if not names:
                    continue
                placeholders = ','.join(['?'] * len(names))
                sql = "SELECT read_name, taxonomy FROM {} WHERE read_name IN ({})".format(
                    READ_NAME_TAXONOMY_TABLE, placeholders)
                for read_name, taxonomy in conn.execute(sql, names).fetchall():
                    to_return[read_name] = [s.strip() for s in taxonomy.split(';')]
        finally:
            conn.close()
        if len(to_return) != len(read_names):
            for r in read_names:
                if r not in to_return:
                    logging.error("Read name {} not found in database.".format(r))
            raise Exception("Found {} read names in metapackage database, expected {}.".format(
                len(to_return), len(read_names)))
        return to_return

    def get_all_taxonomy_strings(self):
        '''Return a list of all taxonomy strings recorded in the database'''
        return [row[0] for row in self._iterate(
            "SELECT DISTINCT taxonomy FROM {}".format(READ_NAME_TAXONOMY_TABLE))]

    def get_taxonomy_by_genome_key(self, genome_keys=None):
        '''Return a dict of genome key (see genome_key_from_read_name, e.g.
        GCF_000744455.1 for a GTDB genome or AMXMAG_0092.fa for a supplemented
        one) to taxonomy string. If genome_keys is given, only those are
        returned.'''
        wanted = set(genome_keys) if genome_keys is not None else None
        key_to_taxonomy = {}
        for read_name, taxonomy in self._iterate(
                "SELECT read_name, taxonomy FROM {}".format(READ_NAME_TAXONOMY_TABLE)):
            key = genome_key_from_read_name(read_name)
            if key is None:
                continue
            if wanted is not None and key not in wanted:
                continue
            if key not in key_to_taxonomy:
                key_to_taxonomy[key] = taxonomy
        return key_to_taxonomy

    def get_marker_counts_of_species(self, taxons):
        '''Return dict of taxonomy string to marker count'''
        taxon_to_count = {}
        conn = self._connect()
        try:
            for taxon_set in iterable_chunks(taxons, 900):
                names = [r for r in taxon_set if r is not None]
                if not names:
                    continue
                placeholders = ','.join(['?'] * len(names))
                sql = "SELECT taxonomy, marker_count FROM {} WHERE taxonomy IN ({})".format(
                    TAXONOMY_MARKER_COUNT_TABLE, placeholders)
                for taxonomy, marker_count in conn.execute(sql, names).fetchall():
                    taxon_to_count[taxonomy] = marker_count
        finally:
            conn.close()
        if len(taxon_to_count) != len(taxons):
            for t in taxons:
                if t not in taxon_to_count:
                    logging.error("Taxon {} not found in database.".format(t))
            raise Exception("Found {} taxons in database, expected {}.".format(
                len(taxon_to_count), len(taxons)))
        return taxon_to_count

    def get_all_marker_counts(self):
        '''Return a dict of all taxonomy strings to marker counts'''
        return {taxonomy: marker_count for taxonomy, marker_count in self._iterate(
            "SELECT taxonomy, marker_count FROM {}".format(TAXONOMY_MARKER_COUNT_TABLE))}

    def _iterate(self, sql, batch_size=50000):
        '''Yield rows from a query in batches (works for both backends), so that
        full-table scans do not buffer the whole result set.'''
        conn = self._connect()
        try:
            cursor = conn.execute(sql)
            while True:
                rows = cursor.fetchmany(batch_size)
                if not rows:
                    break
                for row in rows:
                    yield row
        finally:
            conn.close()
