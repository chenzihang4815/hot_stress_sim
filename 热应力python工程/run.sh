#!/usr/bin/env bash
# Ubuntu / Linux helper. Usage:
#   ./run.sh                         # 3-chip baseline heat+stress
#   ./run.sh cases/case3_4chip.json  # any chips[] JSON
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
if [[ -f .venv/bin/activate ]]; then
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi
JSON="${1:-cases/case3_chips.json}"
python3 scripts/run_case3_stress.py "$JSON"
