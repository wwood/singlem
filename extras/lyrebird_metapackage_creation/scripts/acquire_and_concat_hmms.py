import os
import json
import logging

gpkgs = snakemake.input.gpkgs
output_hmm = snakemake.output.hmm_file
new_hmms_ids = snakemake.output.new_hmms_ids

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

# open CONTENTS.json from each gpkg, extract search_hmms filename, and concatenate to single file
with open(output_hmm, "w+") as w:
    with open(new_hmms_ids, "w+") as w2:
        w2.write("gene\tspkg_name\thmm_filepath\n")
        for gpkg in gpkgs:
            logging.info(f"Reading {gpkg}")
            with open(os.path.join(gpkg, "CONTENTS.json")) as r:
                contents = json.load(r)
                search_hmms = contents.get("search_hmms")
                if search_hmms:
                    with open(os.path.join(gpkg, search_hmms[0])) as r:
                        vog = os.path.basename(gpkg)[:-5]
                        lines = r.readlines()
                        hmm_id = lines[1][:-1].split()[1]
                        for line in lines:
                            w.write(line)
                    w2.write(f"{vog}\t{hmm_id}\t{gpkg}/{search_hmms[0]}\n")
                else:
                    logging.warning(f"No search_hmms found in {gpkg}")

with open(snakemake.output.done, 'w') as _: pass