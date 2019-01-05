import re
import os

class GraftMResult:
    def __init__(self,
                 output_directory,
                 analysing_pairs,
                 search_hmm_files=None):
        self.output_directory = output_directory
        self._search_hmm_files = search_hmm_files
        self._analysing_pairs = analysing_pairs
        if analysing_pairs != False and analysing_pairs != True:
            raise Exception() # To help with refactoring

    def unaligned_sequence_paths(self, require_hits=False):
        '''Return a dict of sample name to hit path, or a Dict of sample name to pair
        of paths if analysing paired data.

        NB: Only return hits when there is at least 1 hit sequence.

        '''
        if not require_hits: raise NotImplementedException()
        paths = {}
        for sample in self.sample_names():
            if self._analysing_pairs:
                paths2 = self.unaligned_paired_sequence_paths_from_sample_name(sample)
                if os.stat(paths2[0]).st_size > 0 or os.stat(paths2[1]).st_size > 0:
                    paths[sample] = paths2
            else:
                path = self.unaligned_sequences_path_from_sample_name(sample)
                if os.stat(path).st_size > 0:
                    paths[sample] = path
        return paths

    def unaligned_sequences_path_from_sample_name(self, sample_name):
        path = os.path.join(self.output_directory,
                            sample_name,
                            "%s_hits.fa" % sample_name)
        return path

    def unaligned_paired_sequence_paths_from_sample_name(self, sample_name):
        paths = [
            os.path.join(self.output_directory,
                         sample_name,
                         'forward',
                         "%s_forward_hits.fa" % sample_name),
            os.path.join(self.output_directory,
                         sample_name,
                         'reverse',
                         "%s_reverse_hits.fa" % sample_name)]
        return paths

    def hmmout_paths_from_sample_name(self, sample_name):
        if self._search_hmm_files is None:
            raise NotImplementedException("Coder needs to grab this from the graftm_package I guess")
        elif len(self._search_hmm_files) == 1:
            return [os.path.join(self.output_directory,
                                 sample_name,
                                 "%s.hmmout.txt" % sample_name)]
        else:
            reg = re.compile(r'\.hmm$')
            hmmout_names = [
                "%s_%s.hmmout.txt" % (reg.sub('', os.path.basename(h)),
                                      sample_name)
                for h in self._search_hmm_files]
            return [os.path.join(
                self.output_directory,
                sample_name,
                hmmout) for hmmout in hmmout_names]

    def sample_names(self, require_hits=False):
        total_list = [f for f in os.listdir(self.output_directory) \
                      if os.path.isdir(os.path.join(self.output_directory, f))]
        if require_hits:
            return [sample for sample in total_list if \
                    os.stat(self.unaligned_sequences_path_from_sample_name(sample)).st_size > 0]
        else:
            return total_list

