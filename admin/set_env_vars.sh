#!/usr/bin/env bash

SINGLEM_METAPACKAGE_PATH_RELATIVE="db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb"
LYREBIRD_METAPACKAGE_PATH_RELATIVE="db/lyrebird_v0.3.0_phrog_v4.1_metapackage_2025_06_23.smpkg.zb"
GTDBTK_DATA_PATH_RELATIVE="db/release226"

# For each of the above variables, check if the file exists in the relative path
# and if it does, set the corresponding environment variable to the absolute path.
for var in SINGLEM_METAPACKAGE_PATH \
           LYREBIRD_METAPACKAGE_PATH \
           GTDBTK_DATA_PATH; do
    relative_path="${var}_RELATIVE"
    abspath="${PIXI_PROJECT_ROOT}/${!relative_path}"
    if [[ -e "${abspath}" ]]; then
        export "$var"="${abspath}"
    else
        echo "File $(basename "${!relative_path}") not found in db/" >&2
    fi
done
