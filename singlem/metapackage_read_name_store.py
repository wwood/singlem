import logging
import tempfile
import extern

from sqlalchemy import create_engine
from sqlalchemy.orm import registry, declarative_base
from sqlalchemy import Column, Integer, String, select

from bird_tool_utils import iterable_chunks

from .singlem_package import SingleMPackage

mapper_registry = registry()
Base = declarative_base()

class MetapackageReadNameStore:
    @staticmethod
    def generate(singlem_package_paths, sqlitedb_path, taxonomy_marker_counts=None):
        engine = create_engine("sqlite+pysqlite:///{}".format(sqlitedb_path),
            echo=logging.getLogger().isEnabledFor(logging.DEBUG),
            future=True)
        ReadNameTaxonomy.metadata.create_all(engine)

        num_packages = 0
        num_read_names = 0

        with tempfile.NamedTemporaryFile(prefix='singlem_metapackage_read_name_store', suffix='.tsv') as temp_file:
            for singlem_package_path in singlem_package_paths:
                singlem_package = SingleMPackage.acquire(singlem_package_path)
                num_packages += 1
                taxonomy_hash = singlem_package.taxonomy_hash()

                for name, taxonomy in taxonomy_hash.items():
                    num_read_names += 1
                    temp_file.write("{}\t{}\t{}\n".format(
                        num_read_names,
                        name,
                        ';'.join(taxonomy)).encode('utf-8'))

            temp_file.flush()

            extern.run("sqlite3 {} '.mode tabs' '.import {} read_name_taxonomy'".format(
                sqlitedb_path, temp_file.name))

            if taxonomy_marker_counts is not None:
                TaxonomyMarkerCount.metadata.create_all(engine)
                with tempfile.NamedTemporaryFile(prefix='singlem_metapackage_read_name_store', suffix='.tsv') as temp_file:
                    for taxonomy, marker_count in taxonomy_marker_counts.items():
                        temp_file.write("{}\t{}\n".format(
                            taxonomy,
                            marker_count).encode('utf-8'))
                    temp_file.flush()

                    extern.run("sqlite3 {} '.mode tabs' '.import {} taxonomy_marker_count'".format(
                        sqlitedb_path, temp_file.name))

        logging.info("Imported {} packages and {} read names.".format(num_packages, num_read_names))

    @staticmethod
    def acquire(sqlitedb_path):
        engine = create_engine("sqlite+pysqlite:///{}".format(sqlitedb_path),
            echo=logging.getLogger().isEnabledFor(logging.DEBUG),
            future=True)

        m = MetapackageReadNameStore()
        m.engine = engine

        return m

    def get_taxonomy_of_reads(self, read_names):
        '''Return dict of read name to taxonomy string'''
        to_return = {}

        for read_name_set in iterable_chunks(read_names, 900): # Must be at least < 1000 for sqlite versions prior to 3.32.0
            names = [r for r in read_name_set if r is not None]
            stmt = select(ReadNameTaxonomy).where(
                ReadNameTaxonomy.read_name.in_(names))
            with self.engine.connect() as conn:
                for res in conn.execute(stmt):
                    to_return[res.read_name] = [s.strip() for s in res.taxonomy.split(';')]
        if len(to_return) != len(read_names):
            raise Exception("Not all read names found in metapackage sqlite3 database")
        return to_return

    def get_all_taxonomy_strings(self):
        '''Return a list of all taxonomy strings recorded in the database'''
        stmt = select(ReadNameTaxonomy.taxonomy).distinct()
        with self.engine.connect() as conn:
            return list([res.taxonomy for res in conn.execute(stmt)])

    def get_marker_counts_of_species(self, taxons):
        '''Return dict of taxonomy string to marker count'''
        taxon_to_count = {}
        for taxon_set in iterable_chunks(taxons, 900): # Must be at least < 1000 for sqlite versions prior to 3.32.0
            names = [r for r in taxon_set if r is not None]
            stmt = select(TaxonomyMarkerCount).where(
                TaxonomyMarkerCount.taxonomy.in_(names))
            with self.engine.connect() as conn:
                for res in conn.execute(stmt):
                    taxon_to_count[res.taxonomy] = res.marker_count
        if len(taxon_to_count) != len(taxons):
            raise Exception("Not all taxons found in metapackage sqlite3 database")
        return taxon_to_count

class ReadNameTaxonomy(Base):
    __tablename__ = "read_name_taxonomy"

    id = Column(Integer, primary_key=True)
    read_name = Column(String, nullable=False, unique=True)
    taxonomy = Column(String, nullable=False)

class TaxonomyMarkerCount(Base):
    __tablename__ = "taxonomy_marker_count"

    id = Column(Integer, primary_key=True)
    taxonomy = Column(String, nullable=False, unique=True)
    marker_count = Column(Integer, nullable=False)