#!/usr/bin/env bash
# Run SVsig-2D (simple model, hg38, FDR 0.1, bin 5e5) locally on every TCGA WGS cohort, from the command line.
#
# Cohorts are discovered from data/input_data/<COHORT>.new.merged.filtered.SV_svsig.csv
# (the files written by data/datasets_to_proper_format_SVsig.ipynb).
# Outputs: results/<COHORT>.new.fdr0.1.bin5e5.hg38.txt  (same naming as the README), logs in results/logs/
#
# Usage:
#   scripts/run_all_cohorts.sh                    # all cohorts, one at a time
#   scripts/run_all_cohorts.sh -j 3               # 3 MATLAB sessions in parallel (watch RAM)
#   scripts/run_all_cohorts.sh BRCA LUAD          # only these cohorts
#   scripts/run_all_cohorts.sh -f                 # force re-run (default skips cohorts with an existing output)
#   scripts/run_all_cohorts.sh -n                 # dry run: list cohorts and commands
# Options: -b BIN_LENGTH (5e5)  -q FDR (0.1)  -g GENOME (hg_38)  -m MATLAB_BIN (matlab)
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
IN_DIR="$REPO/data/input_data"; OUT_DIR="$REPO/results"; LOG_DIR="$OUT_DIR/logs"
BIN=5e5; FDR=0.1; GENOME=hg_38; JOBS=1; FORCE=0; DRY=0; MATLAB="${MATLAB:-matlab}"

while getopts "j:b:q:g:m:fnh" o; do case $o in
  j) JOBS=$OPTARG;; b) BIN=$OPTARG;; q) FDR=$OPTARG;; g) GENOME=$OPTARG;; m) MATLAB=$OPTARG;;
  f) FORCE=1;; n) DRY=1;; h|*) sed -n '2,15p' "$0"; exit 0;; esac; done
shift $((OPTIND-1))

if (( $# )); then COHORTS=("$@"); else
  mapfile -t COHORTS < <(cd "$IN_DIR" && ls *.new.merged.filtered.SV_svsig.csv 2>/dev/null | sed 's/\.new\.merged.*//' | sort)
fi
(( ${#COHORTS[@]} )) || { echo "No cohorts found in $IN_DIR" >&2; exit 1; }
(( DRY )) || command -v "$MATLAB" >/dev/null || { echo "MATLAB not found ('$MATLAB'); use -m /path/to/matlab" >&2; exit 1; }
mkdir -p "$LOG_DIR"
TAG="new.fdr${FDR}.bin${BIN}.${GENOME/_/}"   # e.g. new.fdr0.1.bin5e5.hg38

run_one() {
  local c="$1" sv="$IN_DIR/$1.new.merged.filtered.SV_svsig.csv" out="$OUT_DIR/$1.$TAG.txt"
  [[ -f "$sv" ]] || { echo "[$c] MISSING input $sv"; return 3; }
  if [[ -s "$out" && $FORCE -eq 0 ]]; then echo "[$c] skip (output exists)"; return 0; fi
  local cmd="run2DModel('-wd','$REPO','-sv','$sv','-out','$out','-genome','$GENOME','-bin','$BIN','-fdr','$FDR','-complex','false','-model_exist','false')"
  if (( DRY )); then echo "[$c] $MATLAB -batch \"$cmd\""; return 0; fi
  echo "[$c] start $(date +%T)"
  if "$MATLAB" -batch "$cmd" >"$LOG_DIR/$c.log" 2>&1; then echo "[$c] done  $(date +%T)"
  else echo "[$c] FAILED (see $LOG_DIR/$c.log)"; return 1; fi
}
export -f run_one; export IN_DIR OUT_DIR LOG_DIR BIN FDR GENOME TAG FORCE DRY MATLAB REPO

echo "${#COHORTS[@]} cohorts, $JOBS parallel, bin=$BIN fdr=$FDR genome=$GENOME"
set +e
printf '%s\n' "${COHORTS[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_one {}'
RC=$?; set -e

(( DRY )) && exit 0
echo; echo "== Summary =="
fail=0
for c in "${COHORTS[@]}"; do
  f="$OUT_DIR/$c.$TAG.txt"
  if [[ -s "$f" ]]; then printf '%-12s ok      %6d rows\n' "$c" $(( $(wc -l <"$f") - 1 ))
  else printf '%-12s MISSING/EMPTY\n' "$c"; fail=$((fail+1)); fi
done
echo "$fail of ${#COHORTS[@]} cohorts without output"; exit $(( fail ? 1 : 0 ))
