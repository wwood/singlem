import os
import shutil
import extern

import logging
import pathlib

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

spkg = snakemake.input.singlem_spkg
seqs = snakemake.input.seqs
taxonomy = snakemake.input.taxonomy
off_target_seqs = snakemake.input.off_target_seqs
off_target_taxonomy = snakemake.input.off_target_taxonomy

output_spkg = snakemake.output.spkg

sequence_prefix = snakemake.params.sequence_prefix

if os.path.exists(output_spkg):
    shutil.rmtree(output_spkg)

try:
    print(extern.run(f"singlem regenerate --input-singlem-package {spkg} --output-singlem-package {output_spkg} --min-aligned-percent 0 --input-sequences {seqs} --input-taxonomy {taxonomy} --candidate-decoy-sequences {off_target_seqs} --candidate-decoy-taxonomy {off_target_taxonomy} --sequence-prefix {sequence_prefix} &> {snakemake.log}"))
except Exception as e:
    print(extern.run(f"singlem regenerate --input-singlem-package {spkg} --output-singlem-package {output_spkg} --min-aligned-percent 0 --input-sequences {seqs} --input-taxonomy {taxonomy} --no-further-euks --sequence-prefix {sequence_prefix} &>> {snakemake.log}"))

with open(snakemake.output.done, "w") as _: pass