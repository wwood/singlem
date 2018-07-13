import csv
import logging

class CheckM:
    @staticmethod
    def read_checkm_stats(io, min_completeness, max_contamination):
        # e.g.
        # Bin Id	Marker lineage	# genomes	# markers	# marker sets	0	1	2	3	4	5+	Completeness	Contamination	Strain heterogeneity
        # bin.1	p__Euryarchaeota (UID49)	95	228	153	123	105	0	0	0	0	44.76	0.00	0.00
        next_is_header = True
        ok_bins = set()
        completeness_index = 11
        contamination_index = 12
        for line in csv.reader(io, delimiter="\t"):
            comp = line[completeness_index]
            cont = line[contamination_index]
            if next_is_header:
                if comp != 'Completeness' or cont != 'Contamination':
                    raise Exception("Unexpected checkm csv format detected. "
                                    "This file needs to be generated with CheckM's '--tab_table' flag")
                next_is_header = False
            else:
                completeness = float(comp)
                contamination = float(cont)
                if completeness > 100 or completeness < 0 or contamination < 0:
                    raise Exception("Unexpected entry detected in CheckM line %s" % line)
                if completeness > min_completeness and contamination < max_contamination:
                    logging.debug("Bin %s did pass quality thresholds" % line[0])
                    ok_bins.add(line[0])
                else:
                    logging.debug("Bin %s did not pass quality thresholds" % line[0])
        return ok_bins
