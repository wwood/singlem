#!/usr/bin/env bash
# Activate the pixi dev environment for this repository.

set -e

# Resolve repository root directory
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# pixi uses PIXI_PROJECT_ROOT for activation scripts.
export PIXI_PROJECT_ROOT="$REPO_ROOT"

ENV_PATH="$REPO_ROOT/.pixi/envs/dev"
BIN_PATH="$ENV_PATH/bin"

if [ -d "$BIN_PATH" ]; then
    export PATH="$BIN_PATH:$PATH"
    export CONDA_PREFIX="$ENV_PATH"
else
    echo "Pixi dev environment not found. Run 'pixi install' to create it." >&2
    return 1 2>/dev/null || exit 1
fi

# Run additional activation script to set environment variables
# shellcheck source=/dev/null
. "$REPO_ROOT/admin/set_env_vars.sh"
