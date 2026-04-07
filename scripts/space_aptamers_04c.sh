#!/usr/bin/env bash
# ───────────────────────────────────────────────

set -e  # exit on any error

# ───────────────────────────────────────────────
# Step 1: Create the kmer space file
# ───────────────────────────────────────────────
echo "Generating kmer space file..."
Rscript scripts/kmer_uniques.R 'results' 'graphs'
