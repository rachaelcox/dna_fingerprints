#!/usr/bin/env bash
# space_aptamers_03.sh
# ~~~~~~~~~~~~~~~~~~~~~
# Main pipeline to generate De Bruijn graphs, count k-mers, and combine data.

# Exit on any error
set -e

# ───────────────────────────────
# Step 0 – Input checks
# ───────────────────────────────
if [ $# -lt 1 ]; then
  echo "Usage: $0 <name_of_substrate> [cutoff]"
  echo "Example: $0 streptavidin 10"
  exit 1
fi

name_of_substrate="$1"
cutoff="${2:-20}"           # optional 2nd argument, defaults to 20
num_jobs="${NUM_JOBS:-4}"   # optional env var, defaults to 4 parallel jobs

echo "───────────────────────────────────────────────"
echo "Running pipeline for substrate: ${name_of_substrate}"
echo "Cutoff value: ${cutoff}"
echo "Using ${num_jobs} parallel jobs"
echo "───────────────────────────────────────────────"

# ───────────────────────────────
# Step 1 – De Bruijn graph generation
# ───────────────────────────────
echo " Generating De Bruijn graphs..."
bash ./scripts/05_debruijn_graph_creation.sh "${name_of_substrate}" || {
  echo " Step 1: De Bruijn graphs could not be generated"
  exit 1
}

# ───────────────────────────────
# Step 2 – Count k-mers
# ───────────────────────────────
echo "  Counting k-mers..."
bash ./scripts/06_count_kmers.sh "${name_of_substrate}" || {
  echo " Step 2: K-mers could not be counted"
  exit 1
}

# ───────────────────────────────
# Step 3 – Combine data
# ───────────────────────────────
echo "  Combining data (cutoff = ${cutoff})..."
bash ./scripts/07_combine_data_commands.sh "${name_of_substrate}" "${cutoff}" || {
  echo "  Step 3: Data combination failed"
  exit 1
}

echo -e "  All steps completed successfully for ${name_of_substrate}! (cutoff = ${cutoff}) ♡"
