
It can be useful in some situations to search for sequences in OTU tables. For instance, you may ask "is the most abundant OTU or anything similar in samples B, C or D?" To answer this question make a SingleM database from sample B, C & D's OTU tables:
```
singlem pipe -1 B.fq.gz --otu-table B.otu_table.csv
singlem pipe -1 C.fq.gz --otu-table C.otu_table.csv
singlem pipe -1 D.fq.gz --otu-table D.otu_table.csv
singlem makedb --otu-tables B.otu_table.csv C.otu_table.csv D.otu_table.csv --db BCD.sdb
```
`.sdb` is the conventional file extension for SingleM databases. Then to query this database with windows from sample A:
```
singlem pipe -1 A.fq.gz --otu-table A.otu_table.csv
singlem query --query-otu-table A.otu_table.csv --db BCD.sdb
```
