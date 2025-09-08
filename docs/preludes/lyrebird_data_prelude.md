
The `data` subcommand downloads (or verifies) the reference data used by Lyrebird.
Once it has been downloaded, the environment variable LYREBIRD_METAPACKAGE_PATH
can be used to specify the location of the reference data.

Database versions are documented in the [main Lyrebird documentation page](/Lyrebird).

## Example usage

```
$ lyrebird data --output-directory /path/to/dbs
10/12/2022 01:16:04 PM INFO: Downloading ...
10/12/2022 01:21:49 PM INFO: Extracting files from archive...
10/12/2022 01:22:50 PM INFO: Verifying version and checksums...
10/12/2022 01:23:03 PM INFO: Verification success.
10/12/2022 01:23:03 PM INFO: Finished downloading data
10/12/2022 01:23:03 PM INFO: The environment variable LYREBIRD_METAPACKAGE_PATH can now be set to /path/to/dbs
10/12/2022 01:23:03 PM INFO: For instance, the following can be included in your .bashrc (requires logout and login after inclusion):
10/12/2022 01:23:03 PM INFO: export LYREBIRD_METAPACKAGE_PATH='/path/to/dbs/phrog4.1_v0.2.2024_11_7.smpkg.zb'
```

Then follow the instructions from the final line of above. Set the environment variable so that this data is used automatically, by adding the following to your .bashrc file.

```
export LYREBIRD_METAPACKAGE_PATH='/path/to/dbs/phrog4.1_v0.2.2024_11_7.smpkg.zb'
```
