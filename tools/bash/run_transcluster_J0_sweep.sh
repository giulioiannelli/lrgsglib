#!/usr/bin/env bash
set -euo pipefail

usage(){ cat <<USAGE
Usage: $0 [--start START] [--end END] [--count N] [--backend BACKEND] [--python PYTHON] -- [extra args passed to SCS_TransCluster.py]
Defaults: START=0.0 END=3.0 COUNT=20 BACKEND=cupy PYTHON=python
Example: $0 --start 0 --end 3 --count 20 --backend cupy -- -n a 1000
USAGE
}

START=0.0
END=3.0
COUNT=20
BACKEND="cupy"
PYTHON="python"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2;;
    --end) END="$2"; shift 2;;
    --count) COUNT="$2"; shift 2;;
    --backend) BACKEND="$2"; shift 2;;
    --python) PYTHON="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    --) shift; EXTRA_ARGS+=("$@"); break;;
    *) EXTRA_ARGS+=("$1"); shift;;
  esac
done

if ! command -v "${PYTHON}" &> /dev/null; then
  echo "Warning: python executable '${PYTHON}' not found in PATH."
fi

# Generate J0 values using the selected python executable
mapfile -t J0_VALUES < <("${PYTHON}" - <<PY
import numpy as _np, sys
start = float("${START}")
end = float("${END}")
n = int(${COUNT})
vals = _np.linspace(start, end, n)
sys.stdout.write("\n".join(map(str, vals)))
PY
)

if [[ ${#J0_VALUES[@]} -eq 0 ]]; then
  echo "No J0 values generated. Exiting." >&2
  exit 1
fi

for j0 in "${J0_VALUES[@]}"; do
  echo "--- Running with --J0 ${j0} ---"
  "${PYTHON}" lrgsglib/src/SCS_TransCluster.py 50 --gamma 1.00 --J0 "${j0}" --J 1.00 --g 1.00 --diagonal -1.0 -wd scs_ordParam -sf 20 -na 1000 --backend "${BACKEND}" "${EXTRA_ARGS[@]}"
done

echo "Sweep finished."
