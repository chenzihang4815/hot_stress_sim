#!/usr/bin/env bash
# Ubuntu / Linux runner for Case 2 (die + TIM + spreader + sink).
# Copy this whole hotspot_case2 folder to Linux, then:
#   chmod +x run.sh
#   ./run.sh
#
# HotSpot binary search order:
#   1) $HOTSPOT_BIN
#   2) ./hotspot  (if you copied the compiled binary here)
#   3) ../HotSpot-master/hotspot
#   4) ../HotSpot-master/HotSpot-master/hotspot
#   5) hotspot on $PATH
#
# Compile on Ubuntu if needed:
#   sudo apt-get install -y build-essential
#   cd /path/to/HotSpot-master/HotSpot-master
#   make
#   export HOTSPOT_BIN="$(pwd)/hotspot"

set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")"
NAME="case2_package"
ROWS=64
COLS=64

find_hotspot() {
  local cand
  if [ -n "${HOTSPOT_BIN:-}" ]; then
    printf '%s\n' "$HOTSPOT_BIN"
    return
  fi
  for cand in \
    "./hotspot" \
    "../HotSpot-master/hotspot" \
    "../HotSpot-master/HotSpot-master/hotspot" \
    "../../HotSpot-master/HotSpot-master/hotspot"
  do
    if [ -x "$cand" ]; then
      printf '%s\n' "$cand"
      return
    fi
  done
  if command -v hotspot >/dev/null 2>&1; then
    command -v hotspot
    return
  fi
  return 1
}

if ! HS="$(find_hotspot)"; then
  cat <<'EOF'
[error] hotspot binary not found.

On Ubuntu:
  sudo apt-get install -y build-essential
  cd /path/to/HotSpot-master/HotSpot-master
  make
  export HOTSPOT_BIN="$(pwd)/hotspot"

Then from this folder:
  ./run.sh

Or copy the compiled 'hotspot' binary into this folder and rerun.
EOF
  exit 1
fi

mkdir -p outputs
rm -f outputs/"${NAME}".steady outputs/"${NAME}".grid.steady outputs/"${NAME}".log

echo "[info] hotspot = $HS"
echo "[info] running grid ${ROWS}x${COLS}  (steady only, no transient)"

"$HS" \
  -c "${NAME}.config" \
  -f "${NAME}.flp" \
  -p "${NAME}.ptrace" \
  -materials_file "${NAME}.materials" \
  -model_type grid \
  -grid_rows "$ROWS" \
  -grid_cols "$COLS" \
  -steady_file "outputs/${NAME}.steady" \
  -grid_steady_file "outputs/${NAME}.grid.steady" \
  | tee "outputs/${NAME}.log"

echo
echo "[ok] block temperatures (HotSpot writes Kelvin):  outputs/${NAME}.steady"
if [ -f "outputs/${NAME}.steady" ]; then
  awk '
    BEGIN { print "  name                 T_C" }
    NF>=2 && $1 !~ /^#/ {
      t=$2
      if (t > 200) t = t - 273.15
      printf "  %-18s  %7.3f\n", $1, t
    }
  ' "outputs/${NAME}.steady"
fi
echo "[ok] grid field:  outputs/${NAME}.grid.steady"
echo "Copy outputs/ back to Windows for Python comparison if needed."
