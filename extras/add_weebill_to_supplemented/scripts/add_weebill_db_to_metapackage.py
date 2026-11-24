"""Copy a metapackage, bundle a weebill two-stage database into it, and record
it in CONTENTS.json as a version 7 metapackage. Version 7 loads the version 6
SQLite read-name store unchanged, so nothing else needs to be rebuilt."""

import json
import logging
import os
import shutil

logging.basicConfig(
    filename=snakemake.log[0], level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s')

src = snakemake.input.mpkg
syl2db = snakemake.input.syl2db
dest = snakemake.params.output_metapackage
c = int(snakemake.params.weebill_c)

with open(os.path.join(src, 'CONTENTS.json')) as f:
    contents = json.load(f)
if contents.get('weebill_dbs'):
    raise Exception("Input metapackage already bundles weebill database(s): {}".format(contents['weebill_dbs']))
if contents['singlem_metapackage_version'] < 5:
    raise Exception("Input metapackage version {} is too old".format(contents['singlem_metapackage_version']))

if os.path.exists(dest):
    logging.info("Removing partial output %s", dest)
    shutil.rmtree(dest)
logging.info("Copying %s to %s", src, dest)
# copyfile, not copy2: the destination filesystem may not accept the source's
# ACLs.
shutil.copytree(src, dest, copy_function=shutil.copyfile)

db_name = os.path.basename(syl2db)
logging.info("Copying weebill database %s (c=%d)", syl2db, c)
shutil.copyfile(syl2db, os.path.join(dest, db_name))

contents['singlem_metapackage_version'] = 7
contents['weebill_dbs'] = [{'db': db_name, 'c': c}]
with open(os.path.join(dest, 'CONTENTS.json'), 'w') as f:
    json.dump(contents, f)
logging.info("Wrote CONTENTS.json: %s", {k: v for k, v in contents.items() if k != 'singlem_packages'})

with open(snakemake.output.done, 'w') as f:
    pass
