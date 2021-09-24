import shutil
import logging
import tempfile
import re
import os
import json

from Bio import AlignIO
import extern
from graftm.graftm_package import GraftMPackage

from .singlem_package import SingleMPackage
from .sequence_classes import SeqReader

class PackageHmmTrimmer:
    
    def trim(self, 
        input_package_path,
        output_package_path):

        logging.info("Trimming package {}".format(input_package_path))
        input_package = SingleMPackage.acquire(input_package_path)
        if input_package.window_size() != 60:
            raise Exception("")

        # Copy the entire directory structure to the new output path
        shutil.copytree(input_package_path, output_package_path)
        output_package = SingleMPackage.acquire(output_package_path)

        # Delete the current search HMMs, and make new files, because sometimes
        # the number of search HMMs is different to the number of target domains
        for hmm in output_package.graftm_package().search_hmm_paths():
            os.remove(hmm)

        # Create search HMM for each of Archaea and Bacteria according to target taxonomy
        tax_hash = input_package.taxonomy_hash()
        search_hmm_paths = []
        for target_taxonomy in input_package.target_domains():
            num_hits = 0
            total_hits = 0
            with tempfile.NamedTemporaryFile() as f:
                with open(input_package.graftm_package().unaligned_sequence_database_path()) as seqs:
                    for (name, seq, _) in SeqReader().readfq(seqs):
                        total_hits += 1
                        if tax_hash[name][0].replace('d__','') == target_taxonomy:
                            f.write('>{}\n{}\n'.format(name,seq).encode())
                            num_hits += 1
                logging.info("Found {} of {} hits for domain {}. Generating new search HMM ..".format(num_hits,total_hits,target_taxonomy))
                f.flush()

                # GraftM requires each search HMM basename to be unique, so put
                # the name of the package in the search HMM filename. Also
                # replace dots with underscores as graftm removes everything
                # after the first dot.
                hmm_path = os.path.join(output_package.graftm_package_path(), 
                    "search_{}_{}.hmm".format(
                        target_taxonomy, output_package.graftm_package_basename().replace('.','_')))
                search_hmm_paths.append(os.path.basename(hmm_path))

                extern.run("mafft --thread 8 {} |seqmagick convert --input-format fasta --output-format stockholm - - |hmmbuild --informat stockholm --amino -n {} {} -".format(
                    f.name,
                    "{}.{}".format(input_package.graftm_package_basename(), target_taxonomy),
                    hmm_path
                ))
        # Recreate output graftm package search HMM reference in contents
        output_package.graftm_package()._contents_hash[GraftMPackage.SEARCH_HMM_KEY] = search_hmm_paths
        # save contents file
        with open(output_package.graftm_package().contents_file_path(), 'w') as j:
            json.dump(output_package.graftm_package()._contents_hash, j)

        logging.info("Attempting to create new alignment HMM ..")
        with tempfile.NamedTemporaryFile() as intermediate_hmm:
            # Replace align HMM using all target taxonomy seqs
            example_80char_seq = None
            with tempfile.NamedTemporaryFile(suffix='faa') as align_seqs:
                with open(input_package.graftm_package().unaligned_sequence_database_path()) as f:
                    for (name, seq, _) in SeqReader().readfq(f):
                        if tax_hash[name][0].replace('d__','') in input_package.target_domains():
                            s = '>{}\n{}\n'.format(name,seq)
                            align_seqs.write(s.encode())
                            if example_80char_seq is None and len(seq.replace('-','')) == 80:
                                example_80char_seq = s
                align_seqs.flush()
                extern.run('hmmalign {} {} |hmmbuild --informat stockholm --amino -n {} {} -'.format(
                    input_package.graftm_package().alignment_hmm_path(),
                    align_seqs.name,
                    input_package.graftm_package_basename(),
                    intermediate_hmm.name
                ))

            # Take a sequence from the target taxonomy that is 80aa long, as this is
            # very likely 30aa before, 20aa window, 30aa after. Align the sequence
            # against the align HMM, and extract the position of the 31st and 50th AA.
            with tempfile.NamedTemporaryFile() as f:
                f.write(example_80char_seq.encode())
                f.flush()

                with tempfile.NamedTemporaryFile() as sto:
                    out = extern.run('hmmalign {} {} >{}'.format(
                        intermediate_hmm.name,
                        f.name,
                        sto.name))
                    seq = AlignIO.read(sto.name, 'stockholm')
                    seq = str(seq[0].seq)

                    if example_80char_seq.splitlines()[1][30:50] in seq:
                        logging.info("Appears all good, setting new window")
                        shutil.copyfile(
                            intermediate_hmm.name,
                            output_package.graftm_package().alignment_hmm_path())
                        no_lower_chars = re.sub('[a-z]','',seq)
                        new_window_position = no_lower_chars.index(example_80char_seq.splitlines()[1][30:50])

                        contents_path = output_package.contents_path()
                        output_package._contents_hash[SingleMPackage.SINGLEM_POSITION_KEY] = new_window_position

                        # save contents file
                        with open(os.path.join(output_package.base_directory(), SingleMPackage._CONTENTS_FILE_NAME), 'w') as j:
                            json.dump(output_package._contents_hash, j)
                    else:
                        logging.error("New alignment HMM appears to be different somehow, not replacing alignment HMM")
        
