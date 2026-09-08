
try:
    from importlib.metadata import version as _pkg_version, PackageNotFoundError
    __version__ = _pkg_version("singlem")
except (ImportError, PackageNotFoundError):
    from .version import __version__

OTU_TABLE_OUTPUT_FORMAT = 'standard'
ARCHIVE_TABLE_OUTPUT_FORMAT = 'archive'

CREATE_MIN_ALIGNED_PERCENT = 10
