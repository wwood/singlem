
The `data` subcommand downloads (or verifies) the reference data used by SingleM.
Once it has been downloaded, the environment variable `SINGLEM_METAPACKAGE_PATH`
can be used to specify the location of the reference data.

## Downloading non-default reference data (metpackages)

SingleM currently supports GTDB versions GTDB 07-RS207 onwards, up to the default metapackage version (GTDB R226 as of writing). The newest version of the code is compatible with each of these versions, which correspond to the S4.x and S5.x versions which can be downloaded from [Zenodo](http://dx.doi.org/10.5281/zenodo.5739611) - there each version is a different version of the Zenodo record. 

There also exists a SingleM metapackage for GlobDB, which can be downloaded from the [GlobDB website](https://globdb.org/). As of writing, This metapackage includes many more genomes than that of GTDB.

To use these non-standard reference data, extract the download using `tar -xzf`. The most straightforward way of using the data is probably to use `--metapackage` as a command line argument to `singlem` (rather than specifying the `SINGLEM_METAPACKAGE_PATH` environment variable). The folder to specify as the argument to `--metapackage` is the folder which contains the extracted metapackage folder. The specified folder should either be a SingleM metapackage (canonically ending in `.smpkg`) or a Zenodo Backpack of a SingleM package (canonically ending `.smpkg.zb`).


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
