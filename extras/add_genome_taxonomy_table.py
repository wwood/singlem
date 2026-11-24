#!/usr/bin/env python3
"""Add the genome_taxonomy table to an existing metapackage's read-name store,
so that weebill genomes are mapped to taxonomy by an indexed lookup rather than
a scan of every read name.

Metapackages built by current SingleM already have the table. For one with a
DuckDB store, the table is added in place. For one with a SQLite store (version
3-6, or a version 7 made from one by extras/add_weebill_to_supplemented), a new
DuckDB store with the table is written alongside it and CONTENTS.json is
repointed to it. The SQLite store is then unused and can be deleted.

Modifies the metapackage in place, so run it on a copy if the original must be
kept unchanged.
"""

import argparse
import json
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from singlem.metapackage import Metapackage  # noqa: E402
from singlem.metapackage_read_name_store import MetapackageReadNameStore  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metapackage', required=True, help='Metapackage to modify in place')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s')

    contents_path = os.path.join(args.metapackage, Metapackage._CONTENTS_FILE_NAME)
    with open(contents_path) as f:
        contents = json.load(f)
    if contents[Metapackage.VERSION_KEY] < 3:
        raise Exception("Metapackage version {} has no read-name store".format(contents[Metapackage.VERSION_KEY]))
    store_name = contents[Metapackage.SQLITE_DB_PATH_KEY]
    store_path = os.path.join(args.metapackage, store_name)

    if store_name.endswith('.duckdb'):
        MetapackageReadNameStore.add_genome_taxonomy_table(store_path)
        return

    new_store_name = os.path.splitext(store_name)[0] + '.duckdb'
    new_store_path = os.path.join(args.metapackage, new_store_name)
    logging.info("Converting SQLite store %s to DuckDB store %s", store_path, new_store_path)
    MetapackageReadNameStore.convert_to_duckdb(store_path, new_store_path)
    contents[Metapackage.SQLITE_DB_PATH_KEY] = new_store_name
    tmp_contents_path = contents_path + '.tmp'
    with open(tmp_contents_path, 'w') as f:
        json.dump(contents, f)
    os.replace(tmp_contents_path, contents_path)
    logging.info("Repointed %s to %s. The old store %s is no longer used and can be deleted.",
                 contents_path, new_store_name, store_path)


if __name__ == '__main__':
    main()
