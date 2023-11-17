---
title: SingleM data
---
# singlem data

The `data` subcommand downloads (or verifies) the reference data used by SingleM.
Once it has been downloaded, the environment variable SINGLEM_METAPACKAGE_PATH
can be used to specify the location of the reference data.

## Example usage

```
$ singlem data --output-directory /tmp/dbs
10/12/2022 01:16:04 PM INFO: Downloading ...
10/12/2022 01:21:49 PM INFO: Extracting files from archive...
10/12/2022 01:22:50 PM INFO: Verifying version and checksums...
10/12/2022 01:23:03 PM INFO: Verification success.
10/12/2022 01:23:03 PM INFO: Finished downloading data
10/12/2022 01:23:03 PM INFO: The environment variable SINGLEM_METAPACKAGE_PATH can now be set to /tmp/dbs
10/12/2022 01:23:03 PM INFO: For instance, the following can be included in your .bashrc (requires logout and login after inclusion):
10/12/2022 01:23:03 PM INFO: export SINGLEM_METAPACKAGE_PATH='/tmp/dbs/S3.0.5.metapackage20220806.smpkg.zb/payload_directory'
```

Then follow the instructions from the final line of above. Set the environment variable so that this data is used automatically, by adding the following to your .bashrc file.

```
export SINGLEM_METAPACKAGE_PATH='/tmp/dbs/S3.0.5.metapackage20220806.smpkg.zb/payload_directory'
```

OPTIONS
=======

**\--output-directory** *OUTPUT_DIRECTORY*

  Output directory [required unless SINGLEM_METAPACKAGE_PATH is
    specified]

**\--verify-only**

  Check that the data is up to date and each file has the correct
    checksum

OTHER GENERAL OPTIONS
=====================

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

AUTHORS
=======

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
