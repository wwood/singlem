import logging
from singlem import HmmDatabase

class Chancer:
    def run(self, **kwargs):
        metagenomes = kwargs.pop('metagenomes')
        target_taxonomy = kwargs.pop('target_taxonomy')
        hmmdb = kwargs.pop('hmm_database', HmmDatabase())
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        marker_to_counts = {}
        metagenomes.target_taxonomy = target_taxonomy
        for otu in metagenomes:
            if otu.marker in marker_to_counts:
                marker_to_counts[otu.marker].append(otu.count)
            else:
                marker_to_counts[otu.marker] = [otu.count]

        # Check if every marker gene was detected, warn otherwise
        logging.debug("Found %i different marker genes in the metagenome" %\
                      len(marker_to_counts))

        # Check for unexpected markers
        raise Exception("Not yet implemented")
        import IPython; IPython.embed()



