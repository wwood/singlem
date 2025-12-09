import polars as pl
from csv import DictReader
import logging

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filename=snakemake.log[0])

vcontact_assignments = snakemake.input.vcontact_output_dir + "/exports/final_assignments.csv"
ictv_metadata = snakemake.params.ictv_metadata
metavr_metadata = snakemake.params.metavr_metadata
output_taxonomy = snakemake.output.taxonomy_tsv

genome_to_length = {}

draft_taxa = {}
with open(ictv_metadata, 'r') as ictv_file:
    reader = DictReader(ictv_file, delimiter='\t')
    for row in reader:
        genome = row['genome']
        taxonomy = row['taxonomy']
        draft_taxa[genome] = taxonomy.split(';')
        genome_to_length[genome] = int(row['length'])

with open(metavr_metadata, 'r') as metavr_file:
    reader = DictReader(metavr_file, delimiter='\t')
    for row in reader:
        uvig = row['uvig']
        taxonomy = row['ictv_taxonomy']
        draft_taxa[uvig] = taxonomy.split(';') #TODO: check format
        genome_to_length[uvig] = int(row['length'])

new_taxa = {}
df = pl.DataFrame.read_csv(vcontact_assignments, infer_schema_length=10000)
for genome in draft_taxa.keys():
    draft = draft_taxa[genome]
    vcontact_row = df.filter(pl.col('genome') == genome)
    if vcontact_row.height == 0:
        logging.warning(f"No vContact3 assignment found for genome {genome}, using original taxonomy.")
        new_taxa[genome] = draft
        continue
    tax_levels = ['phylum', 'class', 'order', 'family', 'genus']
    tax_ranks = ['p__', 'c__', 'o__', 'f__', 'g__']
    reference_cols = vcontact_row.select([pl.col(c) for c in df.columns if c.startswith('reference_') and any(lvl in c for lvl in tax_levels)])
    prediction_cols = vcontact_row.select([pl.col(c) for c in df.columns if c.startswith('prediction_') and any(lvl in c for lvl in tax_levels)])

    taxbuilder = []
    for level in tax_levels:
        previous_tax = taxbuilder[-1] if len(taxbuilder) > 0 else None
        ref_col = f'reference_{level}'
        pred_col = f'prediction_{level}'
        ref_tax = reference_cols.select(pl.col(ref_col)).to_series()[0]
        pred_tax = prediction_cols.select(pl.col(pred_col)).to_series()[0]
        if ref_tax is not None and ref_tax != '':
            taxbuilder.append(tax_ranks[tax_levels.index(level)] + ref_tax)
        elif pred_tax is not None and pred_tax != '':
            taxbuilder.append(tax_ranks[tax_levels.index(level)] + pred_tax)
        else:
            if previous_tax is not None:
                last_known = previous_tax.split('__')[1]
            else:
                last_known = 'unclassified'
            taxbuilder.append(tax_ranks[tax_levels.index(level)] + last_known)
    taxbuilder.append('s__' + (draft[6] if len(draft) > 6 and draft[6] != '' else 'unclassified'))
    new_taxa[genome] = taxbuilder

with open(output_taxonomy, 'w') as outfile:
    outfile.write("genome\tlength\ttaxonomy\n")
    for genome, taxonomy in new_taxa.items():
        length = genome_to_length.get(genome, 'NA')
        tax_string = ';'.join(taxonomy)
        outfile.write(f"{genome}\t{length}\t{tax_string}\n")