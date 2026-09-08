
try:
    from importlib.metadata import version as _pkg_version, PackageNotFoundError
    __version__ = _pkg_version("singlem")
except (ImportError, PackageNotFoundError):
    # Either python <3.8, which has no importlib.metadata, or singlem was
    # imported from a source tree that has never been built or installed, so
    # there is no distribution metadata to read. setuptools_scm no longer
    # writes singlem/version.py, so there is nothing to fall back to.
    __version__ = 'unknown'

OTU_TABLE_OUTPUT_FORMAT = 'standard'
ARCHIVE_TABLE_OUTPUT_FORMAT = 'archive'

CREATE_MIN_ALIGNED_PERCENT = 10
