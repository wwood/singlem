import subprocess

class OrfLengthChecker:
    '''GraftM croaks when no ORFs are found that meet the minimum length cutoff.
    This class checks whether input sequence files contain 1 or more ORFs, so we
    don't run into problems later.'''
    @staticmethod
    def check_sequence_file_contains_an_orf(path, min_orf_length):
        '''Returns True or False. Files used here cannot be further used for
        streaming. Only checks the first 1000 lines of the sequence file.'''

        # Cannot use extern here because the SIGPIPE signals generated
        result = subprocess.check_output(['bash','-c',"cat '%s' |zcat --stdout -f  |head -n1000 |orfm -m %i |head -n2" %(
            path, min_orf_length
        )])
        if len(result) == 0:
            return False
        else:
            return True
