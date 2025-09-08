---
title: SingleM data
---
# singlem data

The `data` subcommand downloads (or verifies) the reference data used by SingleM.
Once it has been downloaded, the environment variable `SINGLEM_METAPACKAGE_PATH`
can be used to specify the location of the reference data.

## Downloading non-default reference data (metpackages)

SingleM currently supports GTDB versions GTDB 07-RS207 onwards, up to the default metapackage version (GTDB R226 as of writing). The newest version of the code is compatible with each of these versions, which correspond to the S4.x and S5.x versions which can be downloaded from [Zenodo](http://dx.doi.org/10.5281/zenodo.5739611) - there each version is a different version of the Zenodo record. To use these reference data, extract the download using `tar -xzf` and then either set the environment variable `SINGLEM_METAPACKAGE_PATH` to the path of the extracted directory, or directly use `--metapackage` as a command line argument to `singlem`.


## Example usage

```
$ singlem data --output-directory /path/to/dbs
10/12/2022 01:16:04 PM INFO: Downloading ...
10/12/2022 01:21:49 PM INFO: Extracting files from archive...
10/12/2022 01:22:50 PM INFO: Verifying version and checksums...
10/12/2022 01:23:03 PM INFO: Verification success.
10/12/2022 01:23:03 PM INFO: Finished downloading data
10/12/2022 01:23:03 PM INFO: The environment variable SINGLEM_METAPACKAGE_PATH can now be set to /path/to/dbs
10/12/2022 01:23:03 PM INFO: For instance, the following can be included in your .bashrc (requires logout and login after inclusion):
10/12/2022 01:23:03 PM INFO: export SINGLEM_METAPACKAGE_PATH='/path/to/dbs/S3.0.5.metapackage20220806.smpkg.zb'
```

Then follow the instructions from the final line of above. Set the environment variable so that this data is used automatically, by adding the following to your .bashrc file.

```
export SINGLEM_METAPACKAGE_PATH='/path/to/dbs/S3.0.5.metapackage20220806.smpkg.zb'
```

# OPTIONS

**\--output-directory** *OUTPUT_DIRECTORY*

  Output directory [required unless SINGLEM_METAPACKAGE_PATH is
    specified]

**\--verify-only**

  Check that the data is up to date and each file has the correct
    checksum

# OTHER GENERAL OPTIONS

**\--debug**

  output debug information

**\--version**

  output version information and quit

**\--quiet**

  only output errors

**\--full-help**

  print longer help message

**\--full-help-roff**

  print longer help message in ROFF (manpage) format

# AUTHORS

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
