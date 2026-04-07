#!/bin/bash

# --- SETTINGS ---
# Using the paths you provided
OUTPUT_DIR="blast"
FIG_DIR="figures/images"
QUERY_FASTA="blast/blast_sequences.fasta"
DB_FASTA="consensus/combined_consensus_sequences.fasta"

# Naming BLAST database files
DB_PREFIX="$OUTPUT_DIR/blast_db"
OUTPUT_FILE="$OUTPUT_DIR/blast_results.txt"
IMAGE_PATH="$FIG_DIR/blast_heatmap.png"

# Ensure directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FIG_DIR"

# create blast data base from the combined_consensus_seuqences (contigs)
echo "Indexing database in $OUTPUT_DIR..."
makeblastdb -in "$DB_FASTA" -dbtype nucl -out "$DB_PREFIX" -parse_seqids

# run using the default blast parameters
echo "Starting BLAST search..."
# -outfmt 6 provides a clean tab-separated table
blastn -query "$QUERY_FASTA" \
       -db "$DB_PREFIX" \
       -out "$OUTPUT_FILE" \
       -task blastn \
       -max_target_seqs 1000 \
       -word_size 7 # can also use 11 \
       -reward 2 \
       -penalty -3 \
       -gapopen 5 \
       -gapextend 2 \
       -evalue 0.05 \
       -dust yes \
       -soft_masking true \
       -outfmt 6 # create a clean table \ 
       -num_threads 4

echo "BLAST search complete. Results saved to $OUTPUT_FILE"
echo "Generating heatmap..."

# plot as a heat map using seaborn
python3 << EOF
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap

# Define your custom color palette (Dark Slate to Raspberry Pink)
colors = ["#455a64", "#ad1457"] 
soothing_cmp = LinearSegmentedColormap.from_list("soothing_dark", colors)

cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

try:
    if not os.path.exists("$OUTPUT_FILE") or os.stat("$OUTPUT_FILE").st_size == 0:
        print("No BLAST hits found. Skipping heatmap.")
    else:
        df = pd.read_csv("$OUTPUT_FILE", sep="\t", header=None, names=cols) # read in the blast files
        
        # --- AXES FLIPPED ---
        # index (Rows) = Subjects (Consensus)
        # columns (Cols) = Queries (Baits)
        matrix = df.pivot_table(index='sseqid', columns='qseqid', aggfunc='size', fill_value=0)
        matrix_binary = (matrix > 0).astype(int)

        if matrix_binary.shape[0] >= 1 and matrix_binary.shape[1] >= 1:
            # Clustermap for multiple entries
            if matrix_binary.shape[0] > 1 and matrix_binary.shape[1] > 1:
                g = sns.clustermap(matrix_binary, 
                                   cmap=soothing_cmp, 
                                   metric="jaccard", 
                                   method="complete", 
                                   linewidths=.5, 
                                   linecolor='#263238', 
                                   figsize=(12, 10))
                
                # Format Labels
                plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
                plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
                g.savefig("$IMAGE_PATH")
            
            # Simple Heatmap for single entries
            else:
                plt.figure(figsize=(12, 10))
                sns.heatmap(matrix_binary, cmap=soothing_cmp, cbar=True, linewidths=.5, linecolor='#263238')
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig("$IMAGE_PATH")
            
            print(f"Success! Flipped Heatmap saved to: $IMAGE_PATH")
        else:
            print("Matrix dimensions too small to plot.")

except Exception as e:
    print(f"Python Plotting failed: {e}")
EOF