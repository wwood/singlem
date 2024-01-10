import extern
import re


class HmmLogo:
    @staticmethod
    def hmm_information_content(hmm_path):
        info = extern.run('hmmlogo --no_indel %s' % hmm_path)

        # Parse the output

        # max expected height = 6.45
        # Residue heights
        # 1:  0.037  0.008  0.013  0.019  0.029  0.016  0.009  0.067  0.024  0.119  0.175  0.015  0.010  0.017  0.020  0.026  0.028  0.067  0.005  0.016  ( 0.720)
        # 2:  0.022  0.002  0.025  0.033  0.005  0.013  0.009  0.008  0.052  0.012  0.005  0.054  0.008  0.021  0.025  0.022  0.017  0.011  0.001  0.006  ( 0.353)
        state = 0
        info_content = []
        r = re.compile(r'^ *\d+: .*\( *([0-9.]+)\)$')
        for line in info.split('\n'):
            if line == '':
                continue
            if state == 0:
                if line.startswith('max expected height'):
                    state = 1
                else:
                    raise Exception('Unexpected output from hmmlogo: %s' % line)
            elif state == 1:
                if line.startswith('Residue heights'):
                    state = 2
                else:
                    raise Exception('Unexpected output from hmmlogo: %s' % line)
            elif state == 2:
                if info := r.match(line):
                    info_content.append(float(info.group(1)))
                else:
                    raise Exception('Unexpected output from hmmlogo: %s' % line)
        return info_content
