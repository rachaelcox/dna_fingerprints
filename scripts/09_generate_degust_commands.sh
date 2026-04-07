#!/usr/bin/env bash

# exit on any error
set -e

# ───────────────────────────────────────────────
# Step 0: Parse optional arguments
# ───────────────────────────────────────────────
cutoff=20             # default cutoff value
conditions="0w,1w,3w" # default wash conditions

while [[ "$1" == --* ]]; do
  case "$1" in
    --cutoff)
      if [ -z "$2" ]; then
        echo "Error: --cutoff requires a numeric value"
        exit 1
      fi
      cutoff="$2"
      shift 2
      ;;
    --conditions)
      if [ -z "$2" ]; then
        echo "Error: --conditions requires a comma-separated list"
        exit 1
      fi
      conditions="$2"
      shift 2
      ;;
    *)
      echo "Error: Unknown option '$1'"
      echo "Usage: $0 [--cutoff <value>] [--conditions <cond1,cond2,...>] <exp1> <exp2> ..."
      exit 1
      ;;
  esac
done

# ensure at least one experiment name remains
if [ $# -lt 1 ]; then
  echo "Error: You must specify at least one experiment."
  echo "Usage: $0 [--cutoff <value>] [--conditions <cond1,cond2,...>] <exp1> <exp2> ..."
  exit 1
fi

# ───────────────────────────────────────────────
# Step 1: Setup
# ───────────────────────────────────────────────
mkdir -p cmds logs

for exp in "$@"; do
  mkdir -p results/${exp}/enrichment_nolib
done

# ───────────────────────────────────────────────
# Step 2: Build DEGUST commands
# ───────────────────────────────────────────────
cmd_file="cmds/calculate_enrichment_${cutoff}c_nolib.cmds"
echo "Generating DEGUST commands for cutoff=${cutoff} and conditions=${conditions}..."
: > "${cmd_file}"

for k in {5..15}; do
  for exp in "$@"; do
    counts_file="results/${exp}/counts_combined_nolib/${exp}_${k}mers_${cutoff}c_combined_wide_nolib.csv"
    if [ ! -f "$counts_file" ]; then
      echo "Skipping missing file: ${counts_file}"
      continue
    fi
    echo "Rscript ./scripts/run_degust.R \
  --counts_file ${counts_file} \
  --exp ${conditions} \
  --output_file results/${exp}/enrichment_nolib/${exp}_${k}mers_${cutoff}c_all_washes_enriched_nolib.csv \
  | tee -a logs/${exp}_${k}mers_${cutoff}c_enrichment_nolib.log" \
>> "${cmd_file}"


  done
done

echo "Command list written to ${cmd_file}"

# ───────────────────────────────────────────────
# Step 3: Run in parallel
# ───────────────────────────────────────────────
num_jobs="${NUM_JOBS:-8}"  # default to 8 parallel jobs

echo "🚀 Running DEGUST analysis for experiments: $@"
echo "   → cutoff = ${cutoff}"
echo "   → conditions = ${conditions}"
echo "   → using ${num_jobs} parallel jobs"

cat "${cmd_file}" | parallel --bar --eta -j"${num_jobs}"

echo "All DEGUST enrichment analyses completed successfully ♡"