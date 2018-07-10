import logging
import tempfile

from sequence_extractor import SequenceExtractor
from pipe import SearchPipe

class Renew:
    @staticmethod
    def renew(**kwargs):
        '''Renew an OTU table'''
        input_otus = kwargs.pop('input_otus')
        seqs_list = kwargs['sequences']
        if len(seqs_list) != 1:
            raise Exception("Renew must be run with exactly 1 --sequences argument")
        seqs = seqs_list[0]
        pipe_params = kwargs

        # Collect each of the sequences to extract into a single list
        seq_ids = set()
        for otu in input_otus:
            for name in otu.read_names():
                seq_ids.add(name)
        if len(seq_ids) == 0:
            logging.warn("No reads detected in input OTU tables, not renewing")
            return
        logging.info("Read in {} read names e.g. {}".format(
            len(seq_ids), list(seq_ids)[0]))

        # Extract those sequences
        extracted_seqs = SequenceExtractor().extract_and_read(seq_ids, seqs)

        # Ensure that each and every read was extracted
        if len(seq_ids) != len(extracted_seqs):
            raise Exception("Unexpected number of reads in the OTU table were extracted from the sequence file")

        # Run pipe from the start again on those sequences
        with tempfile.NamedTemporaryFile(suffix='.fna',prefix='singlem_renew') as tf:
            for s in extracted_seqs:
                tf.write(">{}\n{}\n".format(s.name, s.seq))
            tf.flush()

            pipe_params['sequences'] = [tf.name]
            import IPython; IPython.embed()
            SearchPipe().run(**pipe_params)
