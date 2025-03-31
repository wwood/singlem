See top of Snakefile for some instructions on how to run this pipeline.

## Uniprot processing

Download latest Uniprot swissprot annotations

```bash
# Main swissprot annotations - into sprot_to_fasta_and_taxonomy
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

# Running sprot_to_fasta_and_taxonomy.py:
python sprot_to_fasta_and_taxonomy.py -s uniprot_sprot.dat.gz -f uniprot_sprot.fa -t uniprot_sprot_taxonomy.tsv

# Output
# Viruses    17451
# Eukaryota  198986
# Bacteria   336739
# Archaea    19794
```
