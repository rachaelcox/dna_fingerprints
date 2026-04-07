#!/usr/bin/env bash
set -e

# ───────────────────────────────────────────────
# Usage / input check
# ───────────────────────────────────────────────
if [ $# -lt 1 ]; then
  echo "Error: provide at least one substrate/experiment."
  echo "Usage: $0 <substrate1> <substrate2> ..."
  echo "Example: $0 Chitin_Only glass161 streptavidin oligo1"
  exit 1
fi

# All command-line args are substrates
substrates=("$@")

echo "--- Step 1: Combining Data (R) ---"
Rscript scripts/combine_data.R "${substrates[@]}"

echo -e "\n--- Step 2: NMDS and Data Processing (Python) ---"
python3 scripts/nmds_data_processing.py "${substrates[@]}"

echo -e "\n--- Step 3: Figure Creation (R) ---"
Rscript scripts/figures.R

echo -e "\n--- Step 4: Kmer space representations (R) ---"
Rscript scripts/kmer_representation_plot.R 

echo -e "\n--- Figure Generation Complete! ---"
