import re
import os
import csv


class OrfMUtils:
    def un_orfm_name(self, name):
        return re.sub('_\d+_\d+_\d+$', '', name)

class TaxonomyFile:
    def __init__(self, taxonomy_file_path):
        self.sequence_to_taxonomy = {}
        utils = OrfMUtils()
        with open(taxonomy_file_path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.sequence_to_taxonomy[\
                      utils.un_orfm_name(row[0])] = row[1]

    def __getitem__(self, item):
        return self.sequence_to_taxonomy[item]

class HmmDatabase:
    def __init__(self):
        # Array of gpkg names to HmmAndPostion
        self.hmms_and_positions = {}

        for array in [
            ['2.07.ribosomal_protein_L2_rplB.gpkg','DNGNGWU00010_mingle_output_good_seqs.hmm',143],
            ['2.08.ribosomal_protein_L3_rplC.gpkg','DNGNGWU00012_mingle_output_good_seqs.hmm',145],
            ['2.09.ribosomal_protein_L5_rplE.gpkg','DNGNGWU00025_mingle_output_good_seqs.hmm',100],
            ['2.10.ribosomal_protein_L6_rplF.gpkg','DNGNGWU00023_mingle_output_good_seqs.hmm',138],
            ['2.11.ribosomal_protein_L10.gpkg','DNGNGWU00030_mingle_output_good_seqs.hmm',76],
            ['2.12.ribosomal_protein_L11_rplK.gpkg','DNGNGWU00024_mingle_output_good_seqs.hmm',13],
            ['2.13.ribosomal_protein_L14b_L23e_rplN.gpkg','DNGNGWU00014_mingle_output_good_seqs.hmm',82],
            ['2.14.ribosomal_protein_L16_L10E_rplP.gpkg','DNGNGWU00018_mingle_output_good_seqs.hmm',66],
            ['2.15.ribosomal_protein_S2_rpsB.gpkg','DNGNGWU00001_mingle_output_good_seqs.hmm',240],
            ['2.16.ribosomal_protein_S5.gpkg','DNGNGWU00015_mingle_output_good_seqs.hmm',144],
            ['2.17.ribosomal_protein_S7.gpkg','DNGNGWU00017_mingle_output_good_seqs.hmm',84],
            ['2.18.ribosomal_protein_S10_rpsJ.gpkg','DNGNGWU00002_mingle_output_good_seqs.hmm',40],
            ['2.19.ribosomal_protein_S12_S23.gpkg','DNGNGWU00026_mingle_output_good_seqs.hmm',68],
            ['2.20.ribosomal_protein_S15P_S13e.gpkg','DNGNGWU00034_mingle_output_good_seqs.hmm',57],
            ['2.21.ribosomal_protein_S19_rpsS.gpkg','DNGNGWU00016_mingle_output_good_seqs.hmm',33]
          ]:
            hmm_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                            '..', 'db', array[0])
            self.hmms_and_positions[os.path.basename(hmm_path)] = \
                HmmAndPostion(hmm_path,
                               os.path.join(hmm_path, array[1]),
                               array[2]
                               )

    def hmm_paths(self):
        'return an array of absolute paths to the hmms in this database'
        return [hp.hmm_filename for hp in self.hmms_and_positions.values()]

    def gpkg_basenames(self):
        return self.hmms_and_positions.keys()

    def gpkg_paths(self):
        return [h.gpkg_path for _, h in self.hmms_and_positions.iteritems()]

    def __iter__(self):
        for hp in self.hmms_and_positions.values():
            yield hp

class HmmAndPostion:
    def __init__(self, gpkg_path, hmm_filename, best_position):
        self.gpkg_path = gpkg_path
        self.hmm_filename = hmm_filename
        self.best_position = best_position

    def hmm_path(self):
        return os.path.join(self.gpkg_path, self.hmm_filename)
    
    def gpkg_basename(self):
        return os.path.basename(self.gpkg_path)
    
    def hmm_basename(self):
        return os.path.basename(self.hmm_filename)



