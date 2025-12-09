import os
import logging
import argparse
import extern
import polars as pl
import gzip

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download and process MetaVR phages.")
    metadata_exclusive = parser.add_mutually_exclusive_group(required=True)
    metadata_exclusive.description = "Specify either a URL to download the MetaVR metadata or a local file path."
    metadata_exclusive.add_argument("--metadata-download", help="URL to desired MetaVR metadata download. Default: https://www.meta-virome.org/DownloadUvigMetadata")
    metadata_exclusive.add_argument("--metadata-file", type=str, help="Path to local MetaVR metadata file.")
    genome_exclusive = parser.add_mutually_exclusive_group(required=True)
    genome_exclusive.description = "Specify either a URL to download the MetaVR genomes or a local file path."
    genome_exclusive.add_argument("--genome-download", help="URL to desired MetaVR genome download. Default: https://www.meta-virome.org/DownloadUvigSequences")
    genome_exclusive.add_argument("--genome-file", type=str, help="Path to local MetaVR genome file.")
    parser.add_argument("-o", "--outdir", type=str, required=True, help="Output directory for the processed genomes.")
    return parser.parse_args()

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S')

def download_metadata(url: str, outdir: str):
    logging.info(f"Downloading MetaVR metadata from {url}")
    metadata_path = os.path.join(outdir, "metavr_metadata.tsv.gz")
    cmd = f"wget {url} --no-check-certificate -O {metadata_path}"
    extern.run(cmd)
    if not os.path.isfile(metadata_path):
        logging.error("Failed to download MetaVR metadata.")
        return None
    return metadata_path

def download_genomes(url: str, outdir: str):
    logging.info(f"Downloading MetaVR genomes from {url}")
    genomes_path = os.path.join(outdir, "metavr_genomes.fna.gz")
    cmd = f"wget {url} --no-check-certificate -O {genomes_path}"
    extern.run(cmd)
    if not os.path.isfile(genomes_path):
        logging.error("Failed to download MetaVR genomes.")
        return None
    return genomes_path

def filter_metadata(metadata_path) -> pl.DataFrame:
    logging.info("Reading MetaVR metadata into DataFrame")
    df = pl.read_csv(metadata_path, separator="\t", 
                     null_values=["NA", "\\N"], 
                     schema_overrides={"uvig": pl.Utf8, 
                             "length": pl.Int64, 
                             "ictv_taxonomy": pl.Utf8, 
                             "host_taxonomy": pl.Utf8, 
                             "completeness": pl.Float64, 
                             "checkv_contamination": pl.Float64
                            },
                    infer_schema_length=10000)
    # df.select([
    #     pl.col('uvig'),
    #     pl.col('length'),
    #     pl.col('ictv_taxonomy'),
    #     pl.col('host_taxonomy'),
    #     pl.col('completeness'),
    #     pl.col('checkv_contamination'),
    #     pl.col('genome_type')
    # ])
    df = df.filter(pl.col('host_taxonomy').str.contains('Bacteria|Archaea'))
    df = df.filter(pl.col('completeness') >= 80.0)
    df = df.filter(pl.col('checkv_contamination') <= 10)
    df = df.filter(pl.col('genome_type') == 'dsDNA')
    return df

def extract_genomes(genomes_path: str, valid_uvigs: set, outdir: str):
    logging.info("Extracting valid genomes to directory")
    output_genomes_dir = os.path.join(outdir, "genomes")
    os.makedirs(output_genomes_dir, exist_ok=True)
    
    with gzip.open(genomes_path, 'rt') as infile:
        current_uvig_id = None
        current_outfile = None
        
        for line in infile:
            if line.startswith('>'):
                # Close previous file if open
                if current_outfile:
                    current_outfile.close()
                    current_outfile = None
                
                # Extract UVIG ID
                uvig_id = line[1:].split('|')[0]
                current_uvig_id = uvig_id
                
                # Open new file if this UVIG is valid
                if uvig_id in valid_uvigs:
                    output_path = os.path.join(output_genomes_dir, f"{uvig_id}.fna")
                    current_outfile = open(output_path, 'w')
                    current_outfile.write(line)
            else:
                # Write sequence line if we have an open file
                if current_outfile:
                    current_outfile.write(line)
        
        # Close final file if still open
        if current_outfile:
            current_outfile.close()
    
    logging.info(f"Extracted genomes are saved in {output_genomes_dir}")

def main():
    args = parse_arguments()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    if args.metadata_file:
        metadata_path = args.metadata_file
    else:
        metadata_url = args.metadata_download or 'https://www.meta-virome.org/DownloadUvigMetadata'
        metadata_path = download_metadata(metadata_url, outdir)
        if not metadata_path:
            logging.error("Metadata download failed, exiting.")
            return

    if args.genome_file:
        genomes_path = args.genome_file
    else:
        genome_url = args.genome_download or 'https://www.meta-virome.org/DownloadUvigSequences'
        genomes_path = download_genomes(genome_url, outdir)
        if not genomes_path:
            logging.error("Genome download failed, exiting.")
            return

    filtered_df = filter_metadata(metadata_path)
    filtered_df.write_csv(os.path.join(outdir, "metavr_filtered_metadata.tsv"), separator="\t")
    logging.info(f"Filtered metadata saved to {os.path.join(outdir, 'metavr_filtered_metadata.tsv')}")
    valid_uvigs = set(filtered_df['uvig'].to_list())
    extract_genomes(genomes_path, valid_uvigs, outdir)

if __name__ == "__main__":
    main()

