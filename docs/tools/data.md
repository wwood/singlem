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
10/12/2022 01:23:03 PM INFO: export SINGLEM_METAPACKAGE_PATH='/tmp/dbs/S3.0.5.metapackage20220806.smpkg.zb'
```

Then follow the instructions from the final line of above. Set the environment variable so that this data is used automatically, by adding the following to your .bashrc file.

```
export SINGLEM_METAPACKAGE_PATH='/tmp/dbs/S3.0.5.metapackage20220806.smpkg.zb'
```

## Using old GTDB versions

SingleM metapackages based on older versions of GTDB can be downloaded from Zenodo.

| GTDB version | Metapackage version | Compatible SingleM version | Zenodo link |
|--------------|---------------------|----------------------------|-------------|
| R207 | 4.1.0 | 0.18.3 | [https://doi.org/10.5281/zenodo.11107165](https://doi.org/10.5281/zenodo.11107165) |
| R214 | 4.2.2 | 0.18.3 | [https://doi.org/10.5281/zenodo.11123537](https://doi.org/10.5281/zenodo.11123537) |
| R220 | 4.3.0 | 0.18.3 | [https://doi.org/10.5281/zenodo.11323477](https://doi.org/10.5281/zenodo.11323477) |

[zenodo_backpack](https://github.com/centre-for-microbiome-research/zenodo_backpack) can be used to download these files and check their checksums.

```
zenodo_backpack download --doi https://doi.org/10.5281/zenodo.11107165 --output_directory SingleM_metapackage_R207 --bar
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
