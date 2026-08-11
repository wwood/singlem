# Download reference data

Download the reference metapackage used for microbial or phage profiling.

```bash
mkdir -p "$PWD/singlem-data"
singlem data --output-directory "$PWD/singlem-data"
export SINGLEM_METAPACKAGE_PATH="$(find "$PWD/singlem-data" -maxdepth 1 -name '*.smpkg*' -print -quit)"
```

For Lyrebird:

```bash
mkdir -p "$PWD/lyrebird-data"
lyrebird data --output-directory "$PWD/lyrebird-data"
export LYREBIRD_METAPACKAGE_PATH="$(find "$PWD/lyrebird-data" -maxdepth 1 -name '*.smpkg*' -print -quit)"
```
