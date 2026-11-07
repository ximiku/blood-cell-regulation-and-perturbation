#!/bin/sh
set -eu
TASK_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mkdir -p "$TASK_ROOT/.cache/matplotlib" "$TASK_ROOT/.cache/fontconfig"
export MPLCONFIGDIR="$TASK_ROOT/.cache/matplotlib"
export XDG_CACHE_HOME="$TASK_ROOT/.cache"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$TASK_ROOT/.deps"
TASK_PYTHON=${HW2_PYTHON:-/opt/anaconda3/envs/base_clone1/bin/python}
exec "$TASK_PYTHON" "$@"
