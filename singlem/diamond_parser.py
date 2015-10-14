import logging
from graftm.sequence_search_results import DiamondSearchResult, SequenceSearchResult
from singlem import OrfMUtils

class DiamondResultParser:
    def __init__(self, diamond_daa_path):
        self.sequence_to_hit_id = {}
        utils = OrfMUtils()
        logging.debug("Parsing diamond file %s" % diamond_daa_path)
        for arr in DiamondSearchResult.import_from_daa_file(diamond_daa_path).each(\
               [SequenceSearchResult.QUERY_ID_FIELD, SequenceSearchResult.HIT_ID_FIELD]):
            query_id = utils.un_orfm_name(arr[0])
            if query_id in self.sequence_to_hit_id:
                raise Exception("Unexpectedly found duplicate hit in diamond file %s: %s" % (diamond_daa_path, query_id))
            self.sequence_to_hit_id[query_id] = arr[1]
            
        logging.debug("Finished reading diamond file, read in %i assignments" % len(self.sequence_to_hit_id))

    def __getitem__(self, item):
        try:
            return self.sequence_to_hit_id[item]
        except:
            return None #sometimes there is a HMM hit but no DIAMOND hit

