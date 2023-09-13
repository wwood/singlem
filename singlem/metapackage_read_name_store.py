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
    def generate(singlem_package_paths, sqlitedb_path):
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

class ReadNameTaxonomy(Base):
    __tablename__ = "read_name_taxonomy"

    id = Column(Integer, primary_key=True)
    read_name = Column(String, nullable=False, unique=True)
    taxonomy = Column(String, nullable=False)