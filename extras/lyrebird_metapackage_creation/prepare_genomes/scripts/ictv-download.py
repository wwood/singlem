import polars as pl
import argparse
import logging
import os
import extern

### taboos and known typos stolen from https://github.com/christopher-riccardi/ICTVdump
taboo = {'JAEILC010000038', 
        'GECV01031551', 
        'AE006468', 
        'CP015418', 
        'CP000031', 
        'CP000830', 
        'CP001312', 
        'CP001357', 
        'BX897699', 
        'AUXO01792325', 
        'OQ735257', 
        'QQ198719', 
        'QQ198717', 
        'QQ198718',
        'OQ7219011',
        'QNQ73380',
        'QMP84020',
        'QJT73696',
        'QJT73701',
        'QJT73698',
        'SRR2729873', 
        'AUXO017923253', 
        'K03573', 
        'C978956',
        'AF181082',
        'D01221',
        'MT29357',
        'N5326222'} # known invalid accessions

known_typos = {'HQHQ847905':'HQ847905', 
               'EU7257772':'EU725772', 
               'LCM141331':'KX884774', 
               'SKR870013':'KR870013', 
               'AF4362513':'AF362513', 
               'JX4782635':'JX478263',
               'NC010308':'NC_010308',
               'NC021720':'NC_021720',
               'NC035620':'NC_035620',
               'NC028131':'NC_028131',
               'NC027432':'NC_027432',
               'NC012805':'NC_012805',
               'NC012127':'NC_012127'
               } # known typos due to manual annotation, we salvage those by replacing with a list containing the polished accessions

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download and process ICTV phage genomes.")
    parser.add_argument("-u", "--url", default='https://ictv.global/vmr/current', help="URL to desired VMR version. Default: https://ictv.global/vmr/current")
    parser.add_argument("-o", "--outdir", type=str, required=True, help="Output directory for the processed genomes.")
    return parser.parse_args()

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S')

def download_vmr_data(url, outfile):
    logging.info(f"Downloading VMR data from {url}")
    cmd = f"wget {url} --no-check-certificate -O {outfile}"
    extern.run(cmd)
    if not os.path.isfile(outfile):
        logging.error("Failed to download VMR data.")
        raise FileNotFoundError(f"File {outfile} not found after download.")
    logging.info(f"VMR data downloaded to {outfile}")

def read_and_fetch_phage_accessions(vmr_file):
    logging.info("Reading VMR data and fetching phage accessions.")
    xls = pl.read_excel(vmr_file, sheet_id=0)
    vmr_sheet = next((s for s in xls.keys() if 'VMR' in s), None)
    if vmr_sheet is None:
        logging.error("No VMR sheet found in the Excel file.")
        raise ValueError("No VMR sheet found in the Excel file.")
    vmr_df = xls[vmr_sheet]
    logging.info("Reading sheet name: " + vmr_sheet)
    vmr_df = vmr_df.filter(pl.col('Host source').str.contains_any(['bacteria', 'archaea'], ascii_case_insensitive=True))
    vmr_df = vmr_df.filter(pl.col('Genome') == 'dsDNA')
    vmr_df = vmr_df.filter(pl.col('Genome coverage') == 'Complete genome')
    vmr_df = vmr_df.filter(pl.col("Phylum").is_not_null())
    logging.info(f"Filtered to {vmr_df.height} dsDNA phage genomes with complete coverage and known phylum.")
    
    clean_accessions = {}
    for record in vmr_df.iter_rows(named=True):
        acc = record['Virus GENBANK accession'].split()[0]
        if acc in taboo:
            continue
        if acc in known_typos:
            acc = known_typos[acc]
        taxlist = ['Viruses', 
                    record['Phylum'], record['Class'], 
                    record['Order'], record['Family'], 
                    record['Genus'], record['Species'].replace(' ', '_')]
        taxonomy = fix_and_propogate_taxonomy(taxlist)
        clean_accessions[acc] = taxonomy
    logging.info(f"Fetched {len(clean_accessions)} unique phage accessions.")

    return clean_accessions

def fix_and_propogate_taxonomy(taxonomy):
    levels = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    fixed_tax = []
    last_known = 'unclassified'
    for level, name in zip(levels, taxonomy):
        if name is None or len(name.strip()) == 0:
            fixed_tax.append(f"{level}{last_known}")
        else:
            last_known = name
            fixed_tax.append(f"{level}{name}")
    return ';'.join(fixed_tax)

def download_genomes(accessions, output_dir):
    logging.info(f"Preparing to download genomes to {output_dir}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # write accession list to file
    accession_file = os.path.join(output_dir, 'accessions.txt')
    with open(accession_file, 'w') as accfile:
        for acc in accessions.keys():
            accfile.write(f"{acc}\n")
    logging.info(f"Accession list written to {accession_file}")
    
    logging.info("Downloading genomes using efetch.")
    cmd = f"efetch -db nuccore -input {accession_file} -format fasta > {os.path.join(output_dir, 'genomes.fna')}"
    extern.run(cmd)
    logging.info("Genome download completed.")
    
    # split multi-fasta into individual files
    logging.info("Splitting multi-fasta into individual files.")
    genomes_to_taxonomy = {}
    with open(os.path.join(output_dir, 'genomes.fna'), 'r') as infile:
        fasta_data = infile.read().strip().split('>')
    os.makedirs(output_dir + "/genomes", exist_ok=True)
    for entry in fasta_data:
        if entry:
            lines = entry.split('\n')
            header = lines[0].split()[0]
            seq = '\n'.join(lines[1:])
            with open(os.path.join(output_dir, "genomes", f"{header}.fna"), 'w') as outfile:
                outfile.write(f">{header}\n{seq}\n")
            no_version_acc = header.split('.')[0]
            genomes_to_taxonomy[header] = accessions[no_version_acc]
    os.remove(os.path.join(output_dir, 'genomes.fna'))
    logging.info("Individual genome files created.")
    # write taxonomy mapping
    with open(os.path.join(output_dir, 'genome_taxonomy.tsv'), 'w+') as taxfile:
        for genome, taxonomy in genomes_to_taxonomy.items():
            taxfile.write(f"{genome}\t{taxonomy}\n")
    logging.info("Genome taxonomy mapping file created.")

def main():
    args = parse_arguments()
    vmr_file = os.path.join(args.outdir, 'vmr.xlsx')
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    download_vmr_data(args.url, vmr_file)
    accessions = read_and_fetch_phage_accessions(vmr_file)
    download_genomes(accessions, args.outdir)
    logging.info("ICTV phage genome processing completed.")

if __name__ == "__main__":
    main()