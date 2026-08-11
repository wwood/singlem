# Calculate beta diversity

Export OTU sequences and calculate Bray–Curtis distances between samples.

```bash
singlem summarise --input-otu-tables samples.otu.tsv \
  --unifrac-by-otu samples-
convertToEBD.py samples-S3.5.ribosomal_protein_S2_rpsB.unifrac samples.ebd
ExpressBetaDiversity -s samples.ebd -c Bray-Curtis
```
