# Search OTU sequences

Build a database from OTU tables and search it with another sample's OTUs.

```bash
singlem makedb --otu-tables sample-b.otu.tsv sample-c.otu.tsv --db comparison.sdb
singlem query --query-otu-table sample-a.otu.tsv --db comparison.sdb
```
