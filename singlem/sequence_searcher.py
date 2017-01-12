import tempfile
import extern
import string
from uc_file import UCFile
import logging

class SequenceSearcher:
    def global_search(self, query_otu_table_collection,
                     subject_otu_table_collection, cluster_identity):
        '''Search a query OTU table against a subject OTU table, yield over
        UCEntry objects that have been modified so that the query
        and subject are the relevant OtuTableEntry objects rather than
        strings. Or they are None if there are no hits, since
        --output_no_hits is used.

        query_otu_table_collection: OtuTableCollection
        subject_otu_table_collection: OtuTableCollection
        cluster_identity: float or str
            reject hits if have lower identity than this (implemented with vsearch --id).
        '''
        logging.info("Caching query OTUs")
        query_otus = list(query_otu_table_collection)
        logging.info("Caching target OTUs")
        subject_otus = list(subject_otu_table_collection)

        def name_to_index(name):
            return int(string.split(name, ';')[0])

        # write out fasta file numbered to corresponding to the OTU info
        with tempfile.NamedTemporaryFile(prefix='singlem_q_for_vsearch') as query_f:
            for i, u in enumerate(query_otus):
                query_f.write(">%i;size=%i\n" % (i, u.count))
                query_f.write(u.sequence.replace('-','')+"\n")
            query_f.flush()

            with tempfile.NamedTemporaryFile(prefix='singlem_db_for_vsearch') as db_f:
                for i, u in enumerate(subject_otu_table_collection):
                    db_f.write(">%i;size=%i\n" % (i, u.count))
                    db_f.write(u.sequence.replace('-','')+"\n")
                db_f.flush()

                with tempfile.NamedTemporaryFile(prefix='singlem_uc') as uc:
                    command = "vsearch --usearch_global %s --db %s --uc %s --id %s --output_no_hits" % (query_f.name,
                                                                               db_f.name,
                                                                               uc.name,
                                                                               str(cluster_identity))
                    logging.info("Running search")
                    extern.run(command)
                    logging.info("Finished running search")
                    with open(uc.name) as uc_read:
                        for uc_entry in UCFile(uc_read):
                            uc_entry.query = query_otus[name_to_index(uc_entry.query)]
                            if uc_entry.target is not None:
                                uc_entry.target = subject_otus[name_to_index(uc_entry.target)]
                            yield uc_entry
