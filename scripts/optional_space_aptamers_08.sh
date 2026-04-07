#!/usr/bin/env bash
set -euo pipefail

# ───────────────────────────────────────────────
# Step 1: Dendrogram creation and report
# ───────────────────────────────────────────────
echo "Step 1: Clustering and creating dendrogram"
# Output goes to: consensus/clustering/combined_consensus_sequences.fasta.clust.clipkit.clusters
python3 scripts/clustering.py