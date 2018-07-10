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
        output_otu_table = kwargs.pop('otu_table', None)
        archive_otu_table = kwargs.pop('archive_otu_table', None)
        output_extras = kwargs.pop('output_extras')
        singlem_packages = kwargs['singlem_packages']

        pipe_params = kwargs

        # Collect each of the sequences to extract into a single list
        seq_ids = set()
        sample_name = None
        for otu in input_otus:
            if sample_name == None:
                sample_name = otu.sample_name
            elif sample_name != otu.sample_name:
                raise Exception("Renew can currently only be run on OTU tables that contain a single sample")
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
            pipe = SearchPipe()
            new_otus = pipe.run_to_otu_table(**pipe_params)
            sample_name_column = new_otus.fields.index('sample')
            if new_otus is not None:
                for i in range(len(new_otus.data)):
                    new_otus.data[i][sample_name_column] = sample_name
                pipe.write_otu_tables(
                    new_otus,
                    output_otu_table,
                    archive_otu_table,
                    output_extras,
                    singlem_packages)
