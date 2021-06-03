class StreamingHMMSearchResult:
    @staticmethod
    def yield_from_hmmsearch_table(hmmout_path):
        '''yield a query ID fr'''
        # hmmsearch format is
        # qseqid tlen queryname qlen evalue bitscore bias hmmfrom hmmto alifrom alito envfrom envto acc
        #    0    2       3      5     6        7     8      15    16     17     18      19    20   21
        # targetname  accession   tlen queryname accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  hmmfrom   hmmto  alifrom  alito  envfrom envto  acc  description
        #    0           1         2       3        4         5        6       7     8     9  10   11         12       13     14      15       16    17       18       19      20   21      22                                   
        # SRR5161788.1310112_44_5_7 - 86 S2.9.ribosomal_protein_S5.Archaea - 80 8.6e-36 125.9 0.3 1 1 1.2e-38 9.7e-36 125.7 0.3 1 80 5 84 5 84 0.98 -
        # res=HMMSearchResult()
        # res.fields = [
        #                SequenceSearchResult.QUERY_ID_FIELD,
        #                SequenceSearchResult.HMM_NAME_FIELD,
        #                SequenceSearchResult.ACCESSION_ID_FIELD,
        #                SequenceSearchResult.QUERY_LENGTH_FIELD,
        #                SequenceSearchResult.ALIGNMENT_LENGTH_FIELD,
        #                SequenceSearchResult.QUERY_FROM_FIELD,
        #                SequenceSearchResult.QUERY_TO_FIELD,
        #                SequenceSearchResult.HIT_FROM_FIELD,
        #                SequenceSearchResult.HIT_TO_FIELD,
        #                SequenceSearchResult.ALIGNMENT_BIT_SCORE,
        #                SequenceSearchResult.ALIGNMENT_DIRECTION,
        #                ]

        with open(hmmout_path) as f:
            for (i, row) in enumerate([x.rstrip().split() for x in f if not x.startswith('#')]):
                yield row[0]
                # print("row {}, got {}".format(i, ', '.join(row)))
                # alifrom    = int(row[17])
                # alito      = int(row[18])
                # aln_length = (alito-alifrom if alito-alifrom>0 else alifrom-alito)
                # if alito != alifrom: #this actually happens..
                #     res.results.append([row[0],
                #                         row[3],
                #                         row[4],
                #                         row[5],
                #                         aln_length,
                #                         int(row[15]),
                #                         int(row[16]),
                #                         alifrom,
                #                         alito,
                #                         row[7],
                #                         True
                #                         ])