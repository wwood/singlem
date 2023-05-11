from graftm.graftm_package import GraftMPackage, GraftMPackageVersion3
import dendropy
import logging
import tempfile
from Bio import SeqIO
import extern
from .singlem_package import SingleMPackageVersion4
import shutil
import os
import pickle

class PackageCreator:
    def create(self, **kwargs):
        input_graftm_package_path = kwargs.pop('input_graftm_package')
        input_taxonomy = kwargs.pop('input_taxonomy')
        output_singlem_package_path = kwargs.pop('output_singlem_package')
        hmm_position = kwargs.pop('hmm_position')
        window_size = kwargs.pop('window_size')
        target_domains = kwargs.pop('target_domains')
        gene_description = kwargs.pop("gene_description")
        force = kwargs.pop('force')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        if force and os.path.exists(output_singlem_package_path):
            shutil.rmtree(output_singlem_package_path)

        # For protein packages, remove sequences from diamond database that are
        # not in the tree so that hits can be mapped onto the tree and used for
        # alpha and beta diversity metrics.
        gpkg = GraftMPackage.acquire(input_graftm_package_path)
        is_protein_package = SingleMPackageVersion4.graftm_package_is_protein(gpkg)
        logging.info("Detected package type as %s" %
                     ('protein' if is_protein_package else 'nucleotide'))
        if is_protein_package:
            has_tree = True
            try:
                gpkg.reference_package_tree_path()
            except KeyError:
                has_tree = False
                logging.info("Detected tree-less GraftM package.")
            if has_tree:
                tree_leaves = set()
                for node in dendropy.Tree.get(
                        path=gpkg.reference_package_tree_path(),
                        schema='newick').leaf_node_iter():
                    # need to replace here because otherwise they don't line up with the
                    # diamond database IDs
                    node_name = node.taxon.label.replace(' ','_')
                    if node_name in tree_leaves:
                        raise Exception("Found duplicate tree leaf name in graftm package "
                                        "tree. Currently this case is not handled, sorry")
                    tree_leaves.add(node_name)
                for name in tree_leaves: #I don't think there is a 'peek' ?
                    eg_name = name
                    break
                logging.info("Read in %i tree tip names e.g. %s" % (
                    len(tree_leaves), eg_name))

                # Make a new fasta file of all the sequences that are leaves
                found_sequence_names = set()
                num_seqs_unaligned = 0
                filtered_aligned_tempfile = tempfile.NamedTemporaryFile(
                    prefix='singlem_package_creator',
                    suffix='.fasta',
                    mode='w')
                for s in SeqIO.parse(gpkg.unaligned_sequence_database_path(), "fasta"):
                    num_seqs_unaligned += 1
                    if s.id in tree_leaves:
                        if s.id in found_sequence_names:
                            raise Exception("Found duplicate sequence names in graftm unaligned"
                                            " sequence fasta file. Currently this case is not handled,"
                                            " sorry")
                        SeqIO.write([s], filtered_aligned_tempfile, "fasta")
                        found_sequence_names.add(s.id)
                filtered_aligned_tempfile.flush()

                if len(tree_leaves) != len(found_sequence_names):
                    for t in tree_leaves:
                        if t not in found_sequence_names:
                            raise Exception("Found some sequences that were in the tree but not the"
                                            " unaligned sequences database e.g. %s. Something is"
                                            " likely amiss with the input GraftM package" % t)
                    raise Exception("Programming error, shouldn't get here")
                logging.info("All %i sequences found in tree extracted successfully from unaligned"
                            " sequences fasta file, which originally had %i sequences" % (
                                len(found_sequence_names), num_seqs_unaligned))
            else:
                # GraftM package has no tree
                num_seqs_unaligned = 0
                filtered_aligned_tempfile = tempfile.NamedTemporaryFile(
                        prefix='singlem_package_creator',
                        suffix='.fasta',
                        mode='w')
                # okay for now, but copying the file would do just as well
                for s in SeqIO.parse(gpkg.unaligned_sequence_database_path(), "fasta"):
                    num_seqs_unaligned += 1
                    SeqIO.write([s], filtered_aligned_tempfile, "fasta")
                filtered_aligned_tempfile.flush()
                logging.info("All %i sequences successfully copied from unaligned sequences fasta file" % (num_seqs_unaligned))

            # Create a new diamond database
            dmnd_tf = tempfile.NamedTemporaryFile(prefix='singlem_package_creator',suffix='.dmnd')
            cmd = "diamond makedb --in '%s' -d '%s'" % (filtered_aligned_tempfile.name, dmnd_tf.name)
            logging.info("Creating DIAMOND database")
            extern.run(cmd)

        # Create taxonomy hash
        logging.debug("Creating taxonomy hash pickle")
        taxonomy_hash_tf = tempfile.NamedTemporaryFile(prefix='taxonomy_', suffix='.pickle')
        taxonomy_hash_path = taxonomy_hash_tf.name
        tax_hash = {}
        with open(input_taxonomy) as tax:
            for line in tax:
                split_line = line.strip().split('\t')
                tax_hash[split_line[0]] = [taxa.strip() for taxa in split_line[1].split(';')]

        pickle.dump(tax_hash, taxonomy_hash_tf)
        taxonomy_hash_tf.flush()


        # Compile the final graftm/singlem package
        if len(gpkg.search_hmm_paths()) == 1 and \
           gpkg.search_hmm_paths()[0] == gpkg.alignment_hmm_path():
            search_hmms = None
        else:
            search_hmms = gpkg.search_hmm_paths()

        with tempfile.TemporaryDirectory() as tmpdir:
            gpkg_name = os.path.join(
                tmpdir,
                os.path.basename(
                    os.path.abspath(input_graftm_package_path)).replace('.gpkg',''))
            GraftMPackageVersion3.compile(gpkg_name,
                                          gpkg.reference_package_path(),
                                          gpkg.alignment_hmm_path(),
                                          dmnd_tf.name if is_protein_package else None,
                                          gpkg.maximum_range(),
                                          filtered_aligned_tempfile.name if is_protein_package else \
                                              gpkg.unaligned_sequence_database_path(),
                                          gpkg.use_hmm_trusted_cutoff(),
                                          search_hmms)
            logging.debug("Finished creating GraftM package for conversion to SingleM package")
            
            SingleMPackageVersion4.compile(output_singlem_package_path,
                                            gpkg_name, hmm_position, window_size, 
                                            target_domains, gene_description,
                                            taxonomy_hash_path)

            shutil.rmtree(gpkg_name)
            if is_protein_package:
                filtered_aligned_tempfile.close()
                dmnd_tf.close()

            logging.info("SingleM-compatible package creation finished")

        taxonomy_hash_tf.close()
