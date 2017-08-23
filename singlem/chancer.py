import logging
from singlem import HmmDatabase

class Chancer:
    def run_and_print(self, **kwargs):
        print "\t".join(['sample', 'total_seqs', 'homogeneity_index'])
        for sample_prediction in self.predict_samples(**kwargs):
            print sample_prediction

    def predict_samples(self, **kwargs):
        '''Yield a RecoveryPrediction for each sample'''
        metagenomes = kwargs.pop('metagenomes')
        target_taxonomy = kwargs.pop('target_taxonomy') # A list
        hmmdb = kwargs.pop('hmm_database', HmmDatabase())
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        marker_to_counts = {}
        metagenomes.target_taxonomy = target_taxonomy
        last_sample = None
        previous_samples = set()
        for otu in metagenomes:
            # Since we are streaming, check the samples are not intermingled.
            current_sample = otu.sample_name
            if last_sample is None:
                last_sample = current_sample
                previous_samples.add(current_sample)
            elif last_sample != otu.sample_name:
                if current_sample in previous_samples:
                    raise Exception("All OTUs from a particular sample must occur consecutively, found at least one OTU from '%s' that was separate." % current_sample)
                else:
                    previous_samples.add(current_sample)
                    yield self.chance_a_sample(last_sample, hmmdb, marker_to_counts)
                    last_sample = current_sample
                    marker_to_counts = {}

            if otu.marker in marker_to_counts:
                marker_to_counts[otu.marker].append(otu.count)
            else:
                marker_to_counts[otu.marker] = [otu.count]
        if last_sample is not None:
            yield self.chance_a_sample(last_sample, hmmdb, marker_to_counts)

    def chance_a_sample(self, sample_name, hmmdb, marker_to_counts):
        # Check if every marker gene was detected, warn otherwise
        logging.debug("Found %i different marker genes in the metagenome" %\
                      len(marker_to_counts))
        known_markers = set([h.graftm_package_basename() for h in hmmdb.protein_packages()])

        # Check for unexpected markers
        expected_marker_count = 0
        for marker in marker_to_counts.keys():
            if marker not in known_markers:
                logging.info("Marker '%s' is not a current SingleM protein marker, so it is up to you to ensure it is sufficiently single copy" % marker)
            else:
                expected_marker_count += 1
        if expected_marker_count < len(known_markers):
            logging.warn("Found fewer than expected marker genes, suggesting the reported chance of assembly success may be overinflated.")

        per_gene_counts = []
        per_gene_indices = []
        for marker, counts in marker_to_counts.items():
            per_gene_counts.append(sum(counts))
            per_gene_indices.append(max(counts) * Chancer.median(counts) / sum(counts))
        logging.debug("Found elements of prediction: %s and %s" % (
            per_gene_counts, per_gene_indices))

        return RecoveryPrediction(
            sample_name,
            Chancer.mean(per_gene_counts),
            Chancer.mean(per_gene_indices))

    @staticmethod
    def median(numbers):
        s = sorted(numbers)
        n = len(numbers)
        if n % 2 == 0:
            return (float(s[(n-1)/2]) + float(s[n/2])) / 2
        else:
            return s[(n-1)/2]

    @staticmethod
    def mean(numbers):
        return float(sum(numbers)) / max(len(numbers), 1)

class RecoveryPrediction:
    def __init__(self, sample_name, count, homogeneity_index):
        self.sample_name = sample_name
        self.count = count
        self.homogeneity_index = homogeneity_index

    def __str__(self):
        return "\t".join([
            self.sample_name,
            str(self.count),
            str(self.homogeneity_index)])
