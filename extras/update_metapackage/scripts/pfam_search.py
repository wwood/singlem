
# Get around reuqirement for GTDBTK_DATA_PATH
import os
os.environ['GTDBTK_DATA_PATH'] = '/work/microbiome/db/gtdb/gtdb_release207_v2'
# Import the PfamSearch class
from gtdbtk.external.pfam_search import PfamSearch
from gtdbtk.external.pypfam.Scan.PfamScan import PfamScan
from gtdbtk.io.marker.tophit import TopHitPfamFile
import tempfile
import pathlib

pathlib.Path(os.path.dirname(snakemake.output[0])).mkdir(parents=True, exist_ok=True)

pfam_scan = PfamScan(cpu=1, fasta=snakemake.input[0], dir=snakemake.params.pfams)
pfam_scan.search()
pfam_scan.write_results(snakemake.output[0], None, None, None, None)

# with tempfile.NamedTemporaryFile() as tmp:
#     pfam_scan.write_results(tmp.name, None, None, None, None)

#     # Write top hit file
#     tophit_file = TopHitPfamFile(snakemake.output[0], snakemake.wildcards.genome)

#     with open(tmp.name, 'r') as fh_pfam:
#         for line in fh_pfam:
#             if line[0] == '#' or not line.strip():
#                 continue

#             line_split = line.split()
#             gene_id = line_split[0]
#             hmm_id = line_split[5]
#             evalue = float(line_split[12])
#             bitscore = float(line_split[11])
#             tophit_file.add_hit(gene_id, hmm_id, evalue, bitscore)

#     tophit_file.write()

    # # Get top hit
    
    # pfam_search = PfamSearch(threads=1,
    #                         pfam_hmm_dir=snakemake.params.pfams,
    #                         protein_file_suffix=".faa",
    #                         pfam_suffix=".pfam",
    #                         pfam_top_hit_suffix=".tsv",
    #                         checksum_suffix=".sha256",
    #                         output_dir=os.path.dirname(snakemake.output[0]))
    # pfam_search._topHit(tmp.name)

# # List of FASTA files to search
# fasta_files = [[snakemake.wildcards.genome,snakemake.input[0]]]

# # Initialize PfamSearch with desired parameters
# with tempfile.TemporaryDirectory() as tmp_dir:
#     pfam_search = PfamSearch(threads=1,
#                             pfam_hmm_dir=snakemake.params.pfams,
#                             protein_file_suffix=".faa",
#                             pfam_suffix=".pfam",
#                             pfam_top_hit_suffix=".tophits",
#                             checksum_suffix=".sha256",
#                             output_dir=tmp_dir)

#     # Run pfam_search.pl on the list of FASTA files
#     pfam_search.run(fasta_files)

#     import subprocess
#     print(subprocess.check_output(f"ls {tmp_dir}"))
#     shutil.copyfile(f"{tmp_dir}/{snakemake.wildcards.genome}.pfam", "{output}")
