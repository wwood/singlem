#!/usr/bin/env bash

SINGLEM_METAPACKAGE_PATH_RELATIVE="db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb"
LYREBIRD_METAPACKAGE_PATH_RELATIVE="db/phrog4.1_v0.2.2024_11_7.smpkg.zb"

if [[ -e "$SINGLEM_METAPACKAGE_PATH_RELATIVE" ]]; then
    export SINGLEM_METAPACKAGE_PATH="$(realpath "$SINGLEM_METAPACKAGE_PATH_RELATIVE")"
else
    echo "Metapackage $(basename "$SINGLEM_METAPACKAGE_PATH_RELATIVE") not found in db/" >&2
fi

if [[ -e "$LYREBIRD_METAPACKAGE_PATH_RELATIVE" ]]; then
    export LYREBIRD_METAPACKAGE_PATH="$(realpath "$LYREBIRD_METAPACKAGE_PATH_RELATIVE")"
else
    echo "Metapackage $(basename "$LYREBIRD_METAPACKAGE_PATH_RELATIVE") not found in db/" >&2
fi
