#!/usr/bin/env bash
# Run official McPAT 1.0 on an XML description.
# Usage: ./run_mcpat.sh [xml] [print_level]
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
XML="${1:-$HERE/input/ARM_A9_2000.xml}"
LEVEL="${2:-5}"
BIN="$HERE/bin/mcpat"
OUTDIR="$HERE/output"
mkdir -p "$OUTDIR"
base="$(basename "$XML" .xml)"
echo "infile=$XML"
echo "print_level=$LEVEL"
"$BIN" -infile "$XML" -print_level "$LEVEL" -opt_for_clk 1 \
  > "$OUTDIR/${base}.txt" 2> "$OUTDIR/${base}.err"
echo "wrote $OUTDIR/${base}.txt"
python3 "$HERE/scripts/read_mcpat_output.py" "$OUTDIR/${base}.txt"
