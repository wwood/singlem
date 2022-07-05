import os
import json
import logging
import tempfile
import sqlite3
from sqlalchemy import select
import extern

from singlem.sequence_database import SequenceDatabase
from singlem_database_models import *

class TaxonomyDB:
    version = 1
    SQLITE_DB_NAME = 'taxonomy.sqlite3'

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    VERSION_KEY = 'singlem_taxonomy_database_version'

    _REQUIRED_KEYS = {1: [VERSION_KEY]}

    @staticmethod
    def create_from_otu_table(input_db, output_db):
        # Acquire sdb
        logging.info("Acquiring sequence database ..")
        sdb = SequenceDatabase.acquire(input_db)

        # Make new taxdb folder and CONTENTS.json
        # ensure db does not already exist
        if os.path.exists(output_db):
            raise Exception("Cowardly refusing to overwrite already-existing taxonomy database path '%s'" % output_db)
        logging.info("Creating SingleM taxonomy database at {}".format(output_db))
        os.makedirs(output_db)

        # Create contents file
        contents_file_path = os.path.join(output_db, TaxonomyDB._CONTENTS_FILE_NAME)
        with open(contents_file_path, 'w') as f:
            json.dump({
                TaxonomyDB.VERSION_KEY: 1,
            }, f)
        
        def lca_taxonomy(taxonomy_strings):
            lca = []
            hit_taxonomies = list([list(t.split('; ')) for t in taxonomy_strings])
            for (i, taxon) in enumerate(hit_taxonomies[0]):
                if all([len(h) > i and h[i]==taxon for h in hit_taxonomies]):
                    lca.append(taxon)
                else:
                    break
            if lca == []:
                return 'Root'
            else:
                return '; '.join(lca)

        # For each OTU sequence, assign it taxonomy by taking the LCA of the associated sequences
        # Write output to tsv
        with tempfile.NamedTemporaryFile() as f:
            # Query for nucleotide sequences including taxonomy column of otus table
            with sdb.engine.connect() as conn:
                batch_size = 10000
                builder = select(Otu.id, Otu.marker_id, Otu.marker_wise_sequence_id, Otu.taxonomy_id
                    ).order_by(Otu.marker_id, Otu.marker_wise_sequence_id
                    ).execution_options(yield_per=batch_size)
                last_sequence_id = None
                last_marker_id = None
                prev_taxonomies = []
                count = 0
                for batch in conn.execute(builder).partitions(batch_size):
                    for row in batch:
                        if last_sequence_id == None:
                            last_sequence_id = row.marker_wise_sequence_id
                            last_marker_id = row.marker_id
                            prev_taxonomies = [sdb.get_taxonomy_via_cache(row.taxonomy_id)]
                        elif last_marker_id != row.marker_id or \
                            last_sequence_id != row.marker_wise_sequence_id:
                            f.write(("\t".join((
                                str(count+1),
                                str(last_marker_id),
                                str(last_sequence_id),
                                lca_taxonomy(prev_taxonomies)
                            ))+"\n").encode())
                            count += 1
                            prev_taxonomies = [sdb.get_taxonomy_via_cache(row.taxonomy_id)]
                            last_sequence_id = row.marker_wise_sequence_id
                            last_marker_id = row.marker_id
                        else:
                            prev_taxonomies.append(sdb.get_taxonomy_via_cache(row.taxonomy_id))
                    # Process the last one
                if last_sequence_id is not None:
                    f.write(("\t".join((
                        str(count+1),
                        str(last_marker_id),
                        str(last_sequence_id),
                        lca_taxonomy(prev_taxonomies)
                    ))+"\n").encode())
                    count += 1
                logging.info("Collected info on {} sequences".format(count))
                
                f.flush()

                logging.info("Creating taxon query table and importing data ..")
                sqlite_db_path = os.path.join(output_db, TaxonomyDB.SQLITE_DB_NAME)
                extern.run('sqlite3 {}'.format(sqlite_db_path), stdin= \
                    "CREATE TABLE window_taxonomy (id INTEGER PRIMARY KEY AUTOINCREMENT,"
                    " marker_id INTEGER,"
                    " marker_wise_sequence_id INTEGER,"
                    " taxonomy text);\n"
                    '.separator "\\t"\n'
                    ".import {} window_taxonomy".format(f.name))
                logging.info("Finished import")

                logging.info("Creating SQL indexes ..")
                logging.debug("Connecting to db %s" % sqlite_db_path)
                db = sqlite3.connect(sqlite_db_path)
                c = db.cursor()
                c.execute("CREATE INDEX window_taxonomy_marker_id on window_taxonomy (marker_id)")
                c.execute("CREATE INDEX window_taxonomy_marker_wise_sequence_id on window_taxonomy (marker_wise_sequence_id)")
                c.execute("CREATE UNIQUE INDEX window_taxonomy_marker_and_sequence_id on window_taxonomy (marker_id, marker_wise_sequence_id)")
                db.commit()
                logging.info("Finished")
                        

