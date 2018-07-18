import json
import os
import logging
from graftm.graftm_package import GraftMPackage
import shutil
import hashlib

class InsufficientSingleMPackageException(Exception): pass
class MalformedSingleMPackageException(Exception): pass

class SingleMPackage:
    '''Class to represent a SingleM package. To start using a package, run

    pkg = SingleMPackage.acquire('/path/to/a.spkg')

    and then retrieve values with e.g.

    pkg.graftm_package_path()
      #=> '/path/to/a.spkg/the.gpkg'
    '''

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    # The key names are unlikely to change across package format versions,
    # so store them here in the superclass
    GRAFTM_PACKAGE_KEY = 'graftm_package_path'
    VERSION_KEY = 'singlem_package_version'
    SINGLEM_POSITION_KEY = 'singlem_hmm_position'
    SINGLEM_WINDOW_SIZE_KEY = 'singlem_window_size'
    ALIGNMENT_HMM_SHA256_KEY = 'alignment_hmm_sha256'
    SINGLEM_PACKAGE_SHA256_KEY = 'singlem_package_sha256'

    _CURRENT_FORMAT_VERSION = 1

    _REQUIRED_KEYS = {'1': [
                             VERSION_KEY,
                             GRAFTM_PACKAGE_KEY,
                             SINGLEM_POSITION_KEY,
                             ALIGNMENT_HMM_SHA256_KEY
                             ],
                      '2': [
                             VERSION_KEY,
                             GRAFTM_PACKAGE_KEY,
                             SINGLEM_POSITION_KEY,
                             SINGLEM_WINDOW_SIZE_KEY,
                             ALIGNMENT_HMM_SHA256_KEY
                             ],
                      }


    @staticmethod
    def acquire(singlem_package_path):
        '''Acquire a new singlem Package

        Parameters
        ----------
        singlem_package_path: str
            path to base directory of singlem package
        '''

        with open(os.path.join(
                singlem_package_path,
                SingleMPackage._CONTENTS_FILE_NAME)) as f:
            contents_hash = json.load(f)

        v=contents_hash[SingleMPackage.VERSION_KEY]
        logging.debug("Loading version %i SingleM package: %s" % (v, singlem_package_path))
        if v == 1:
            pkg = SingleMPackageVersion1()
        elif v == 2:
            pkg = SingleMPackageVersion2()
        else:
            raise InsufficientSingleMPackageException("Bad SingleM package version: %s" % str(v))

        pkg._contents_hash = contents_hash
        pkg._base_directory = singlem_package_path
        # check we are at current version otherwise choke
        pkg.check_universal_keys(v)
        pkg.check_required_keys(SingleMPackage._REQUIRED_KEYS[str(v)])
        return pkg

    def check_universal_keys(self, version):
        h = self._contents_hash
        try:
            v = h[SingleMPackage.VERSION_KEY]
        except KeyError:
            raise MalformedSingleMPackageException("No version information in graftm package")
        if v != version:
            raise MalformedSingleMPackageException("Bad version: %s" % v)

    def check_required_keys(self, required_keys):
        '''raise InsufficientGraftMPackageException if this package does not
        conform to the standard of the given package'''
        h = self._contents_hash
        for key in required_keys:
            if key not in h:
                raise MalformedSingleMPackageException("Package missing key %s" % key)

    def __getitem__(self, key):
        '''Return the value of the given key from the contents file'''
        return self.contents_hash[key]

    def contents_path(self):
        return os.path.join(self._base_directory, SingleMPackage._CONTENTS_FILE_NAME)

    def base_directory(self):
        return self._base_directory

class SingleMPackageVersion1(SingleMPackage):
    version = 1 # don't change me bro

    def __init__(self):
        self.graftm_package_cache = None

    def graftm_package_path(self):
        return os.path.join(self._base_directory,
                            self._contents_hash[SingleMPackage.GRAFTM_PACKAGE_KEY])

    def graftm_package(self):
        if self.graftm_package_cache is None:
            self.graftm_package_cache = GraftMPackage.acquire(self.graftm_package_path())
        return self.graftm_package_cache

    def graftm_package_basename(self):
        return os.path.basename(self.graftm_package_path())

    def singlem_position(self):
        return self._contents_hash[SingleMPackage.SINGLEM_POSITION_KEY]

    def alignment_hmm_sha256(self):
        return self._contents_hash[SingleMPackage.ALIGNMENT_HMM_SHA256_KEY]

    def singlem_package_sha256(self):
        return self._contents_hash[SingleMPackage.SINGLEM_PACKAGE_SHA256_KEY]

    def hmm_path(self):
        return self.graftm_package().alignment_hmm_path()

    def hmm_basename(self):
        return os.path.basename(self.hmm_path())

    def calculate_alignment_hmm_sha256(self):
        return hashlib.sha256(open(self.graftm_package().alignment_hmm_path()).read()).hexdigest()

    def calculate_singlem_package_sha256(self):
        h = hashlib.sha256()
        # Maybe the reference package contents should be cached instead
        # of the reference package contents file so that randomly
        # generated files don't interfere, but eh for now.
        h.update(str(self.graftm_package()._refpkg_contents()))
        h.update(str(self.version))
        h.update(str(self.singlem_position()))
        files_to_hash = [self.graftm_package().alignment_hmm_path()]
        if self.is_protein_package():
            files_to_hash.append(self.graftm_package().diamond_database_path())
        files_to_hash.append(self.graftm_package().unaligned_sequence_database_path())
        files_to_hash += self.graftm_package().search_hmm_paths()
        for f in files_to_hash:
            h.update(open(f).read())
        return h.hexdigest()

    def window_size(self):
        return 60

    @staticmethod
    def graftm_package_is_protein(graftm_package):
        return graftm_package.is_protein_package()

    def is_protein_package(self):
        '''Return true if this package is an Amino Acid alignment package, otherwise
        False i.e. it is a nucleotide package.

        '''
        if not hasattr(self, '_is_protein_package'):
            self._is_protein_package = SingleMPackageVersion1.graftm_package_is_protein(
                self.graftm_package())
        return self._is_protein_package

    @staticmethod
    def compile(output_package_path, graftm_package_path, singlem_position):
        '''Create a new SingleM package with the given inputs. Any files
        specified as parameters are copied into the final package so can
        be removed after calling this function.

        Parameters
        ----------
        output_package_path: str
            path to the package being created (must not exist)
        graftm_package_path: str
            path to graftm package internal to the singlem package
        singlem_position: int
            the position in the HMM where the SingleM window starts

        Returns
        -------
        Nothing
        '''

        if os.path.exists(output_package_path):
            raise Exception("Not writing new SingleM package to already existing file/directory with name %s" % output_package_path)
        os.mkdir(output_package_path)

        graftm_package = GraftMPackage.acquire(graftm_package_path)
        if graftm_package.version != 3:
            raise Exception("SingleM packages can only be created from version 3 GraftM packages at this point.")
        graftm_package_basename = os.path.basename(
            output_package_path.replace('.spkg','').replace('.gpkg',''))
        logging.info("Using GraftM package name %s" % graftm_package_basename)
        if graftm_package_basename == SingleMPackage._CONTENTS_FILE_NAME:
            raise Exception("Name of GraftM package cannot be %s" % SingleMPackage._CONTENTS_FILE_NAME)
        shutil.copytree(graftm_package_path, os.path.join(output_package_path, graftm_package_basename))

        singlem_package = SingleMPackageVersion1()
        singlem_package._contents_hash = {SingleMPackage.VERSION_KEY: singlem_package.version,
                                          SingleMPackage.GRAFTM_PACKAGE_KEY: graftm_package_basename,
                                          SingleMPackage.SINGLEM_POSITION_KEY: singlem_position
                                          }
        singlem_package._base_directory = output_package_path

        # calculate the sha256 values
        singlem_package._contents_hash[SingleMPackage.ALIGNMENT_HMM_SHA256_KEY] = \
            singlem_package.calculate_alignment_hmm_sha256()
        singlem_package._contents_hash[SingleMPackage.SINGLEM_PACKAGE_SHA256_KEY] = \
            singlem_package.calculate_singlem_package_sha256()

        # save contents file
        json.dump(singlem_package._contents_hash,
                  open(os.path.join(output_package_path, SingleMPackage._CONTENTS_FILE_NAME), 'w'))

class SingleMPackageVersion2(SingleMPackageVersion1):
    '''Version 2 packages are like version 1 packages except that the window size
    is defined.'''
    version = 2 # don't change me bro

    def window_size(self):
        return self._contents_hash[SingleMPackage.SINGLEM_WINDOW_SIZE_KEY]

    @staticmethod
    def compile(output_package_path, graftm_package_path, singlem_position, window_size):
        if os.path.exists(output_package_path):
            raise Exception("Not writing new SingleM package to already existing file/directory with name %s" % output_package_path)
        os.mkdir(output_package_path)

        graftm_package = GraftMPackage.acquire(graftm_package_path)
        if graftm_package.version != 3:
            raise Exception("SingleM packages can only be created from version 3 GraftM packages at this point.")
        # Use abspath before basename so that trailing slashes are dealt with.
        graftm_package_basename = os.path.basename(
            os.path.abspath(output_package_path).replace('.spkg','').replace('.gpkg',''))
        logging.info("Using GraftM package name %s" % graftm_package_basename)
        if graftm_package_basename == SingleMPackage._CONTENTS_FILE_NAME:
            raise Exception("Name of GraftM package cannot be %s" % SingleMPackage._CONTENTS_FILE_NAME)
        shutil.copytree(graftm_package_path, os.path.join(output_package_path, graftm_package_basename))

        singlem_package = SingleMPackageVersion2()
        singlem_package._contents_hash = {SingleMPackage.VERSION_KEY: singlem_package.version,
                                          SingleMPackage.GRAFTM_PACKAGE_KEY: graftm_package_basename,
                                          SingleMPackage.SINGLEM_POSITION_KEY: singlem_position,
                                          SingleMPackage.SINGLEM_WINDOW_SIZE_KEY: window_size
                                          }
        singlem_package._base_directory = output_package_path

        if singlem_package.is_protein_package() and window_size % 3 != 0:
            raise Exception("For protein packages, the window size must be specified in base pairs. However, the window_size specified is not divisible by 3.")

        # calculate the sha256 values
        singlem_package._contents_hash[SingleMPackage.ALIGNMENT_HMM_SHA256_KEY] = \
            singlem_package.calculate_alignment_hmm_sha256()
        singlem_package._contents_hash[SingleMPackage.SINGLEM_PACKAGE_SHA256_KEY] = \
            singlem_package.calculate_singlem_package_sha256()

        # save contents file
        with open(os.path.join(
                output_package_path, SingleMPackage._CONTENTS_FILE_NAME), 'w') as f:
            json.dump(singlem_package._contents_hash, f)
