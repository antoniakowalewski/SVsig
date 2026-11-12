#!/usr/bin/env bash
# Run SVsig-2D (simple model, hg38, FDR 0.1, bin 5e5) locally on every TCGA WGS cohort, from the command line.
#
# Step 1: convert raw <COHORT>.new.merged.filtered.SV.tsv files (-i RAW_DIR) to SVsig CSVs with
#         data/datasets_to_proper_format_SVsig.py (skipped with -s, or if RAW_DIR does not exist
#         but data/input_data already has converted files).
# Step 2: run SVsig on every cohort discovered in data/input_data/<COHORT>.new.merged.filtered.SV_svsig.csv
# Outputs: results/<COHORT>.new.fdr0.1.bin5e5.hg38.txt  (same naming as the README), logs in results/logs/
#
# Usage:
#   scripts/run_all_cohorts.sh                    # all cohorts, one at a time
#   scripts/run_all_cohorts.sh -j 3               # 3 MATLAB sessions in parallel (watch RAM)
#   scripts/run_all_cohorts.sh BRCA LUAD          # only these cohorts
#   scripts/run_all_cohorts.sh -f                 # force re-run (default skips cohorts with an existing output)
#   scripts/run_all_cohorts.sh -n                 # dry run: list cohorts and commands
#   scripts/run_all_cohorts.sh -i /path/to/NEW_filtered   # raw SV .tsv folder for the conversion step
#   scripts/run_all_cohorts.sh -s                 # skip the conversion step
# Options: -i RAW_DIR -p PYTHON(python3) -b BIN_LENGTH (5e5)  -q FDR (0.1)  -g GENOME (hg_38)  -m MATLAB_BIN (matlab)
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
IN_DIR="$REPO/data/input_data"; OUT_DIR="$REPO/results"; LOG_DIR="$OUT_DIR/logs"
BIN=5e5; FDR=0.1; GENOME=hg_38; JOBS=1; RAW_DIR="${SVSIG_RAW_DIR:-}"; SKIP_CONVERT=0; PYTHON="${PYTHON:-python3}"; FORCE=0; DRY=0; MATLAB="${MATLAB:-matlab}"

while getopts "j:b:q:g:m:i:p:sfnh" o; do case $o in
  j) JOBS=$OPTARG;; b) BIN=$OPTARG;; q) FDR=$OPTARG;; g) GENOME=$OPTARG;; m) MATLAB=$OPTARG;; i) RAW_DIR=$OPTARG;; p) PYTHON=$OPTARG;; s) SKIP_CONVERT=1;;
  f) FORCE=1;; n) DRY=1;; h|*) sed -n '2,20p' "$0"; exit 0;; esac; done
shift $((OPTIND-1))

# ---- Step 1: raw SV tables -> SVsig input CSVs ----
CONVERT="$REPO/data/datasets_to_proper_format_SVsig.py"
if (( ! SKIP_CONVERT )); then
  conv_args=(-o "$IN_DIR"); [[ -n "$RAW_DIR" ]] && conv_args+=(-i "$RAW_DIR")
  if (( DRY )); then echo "[convert] $PYTHON $CONVERT ${conv_args[*]} $*"
  else
    set +e; "$PYTHON" "$CONVERT" "${conv_args[@]}" "$@"; rc=$?; set -e
    if (( rc == 2 )) && compgen -G "$IN_DIR/*.new.merged.filtered.SV_svsig.csv" >/dev/null; then
      echo "[convert] raw SV folder not found; using existing converted files in $IN_DIR" >&2
    elif (( rc != 0 )); then echo "[convert] FAILED (exit $rc). Use -i RAW_DIR, or -s to skip." >&2; exit 1; fi
  fi
fi

if (( $# )); then COHORTS=("$@"); else
  mapfile -t COHORTS < <(cd "$IN_DIR" 2>/dev/null && ls *.new.merged.filtered.SV_svsig.csv 2>/dev/null | sed 's/\.new\.merged.*//' | sort)
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
