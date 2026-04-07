#!/usr/bin/env bash
# ───────────────────────────────────────────────

set -e  # exit on any error

# ───────────────────────────────────────────────
# Step 0: Parse and forward arguments
# ───────────────────────────────────────────────

args=("$@")

if [ ${#args[@]} -lt 1 ]; then
  echo "Error: You must specify at least one experiment."
  echo "Usage: $0 [--cutoff <value>] [--conditions <cond1,cond2,...>] <exp1> <exp2> ..."
  exit 1
fi

# ───────────────────────────────────────────────
# Step 2: Run library subtraction stage
# ───────────────────────────────────────────────
echo "───────────────────────────────────────────────"
echo "Step 1: Subtract library kmers (08_subtract_libraries.sh)"
echo "───────────────────────────────────────────────"

bash "./scripts/08_subtract_libraries.sh" "${args[@]}"

echo "Library subtraction completed successfully"
echo ""

# ───────────────────────────────────────────────
# Step 3: Run DEGUST enrichment analysis stage
# ───────────────────────────────────────────────
echo "───────────────────────────────────────────────"
echo "Step 2: Run DEGUST enrichment analysis (09_generate_degust_commands.sh)"
echo "───────────────────────────────────────────────"

bash "./scripts/09_generate_degust_commands.sh" "${args[@]}"

echo ""
echo "Full pipeline completed successfully ♡"