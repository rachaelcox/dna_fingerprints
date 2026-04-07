#!/usr/bin/env bash
set -euo pipefail

# ───────────────────────────────────────────────
# Step 0: Input check
# ───────────────────────────────────────────────
if [ $# -lt 1 ]; then
  echo "Error: you must specify at least one experiment (e.g. oligo1 oligo2 ...)"
  echo "Usage: $0 [NUM_JOBS] <exp1> <exp2> ..."
  exit 1
fi

# ───────────────────────────────────────────────
# Step 0.5: NUM_JOBS (default = 8)
# ───────────────────────────────────────────────
NUM_JOBS=8

# if first argument is an integer, treat it as NUM_JOBS
if [[ "$1" =~ ^[0-9]+$ ]]; then
  NUM_JOBS="$1"
  shift
fi

# ───────────────────────────────────────────────
# Step 1: Setup output dirs
# ───────────────────────────────────────────────
mkdir -p cmds

for exp in "$@"; do
  mkdir -p "results/${exp}/consensus_seqs_nolib"
done

# # ───────────────────────────────────────────────
# # Step 2: Build command file
# # ───────────────────────────────────────────────
 cmd_file="cmds/top_kmer_weights_nolib.cmds"
 : > "${cmd_file}"

 for exp in "$@"; do
   for k in {5..15}; do
     echo "python3 ./scripts/calc_consensus_seq.py \
 --input_file results/${exp}/consensus_seqs_nolib/${exp}_${k}mers_weights.csv \
 --calc_n 100 \
 --edges" >> "${cmd_file}"
   done
 done

 echo "Wrote commands to: ${cmd_file}"

# ───────────────────────────────────────────────
# Step 3: Run
# ───────────────────────────────────────────────
echo "Running with GNU parallel (-j${NUM_JOBS})..."

# --joblog: creates a tab-separated file with exit codes and commands
# --resume: allows you to re-run the script and it will only pick up failed jobs
parallel -j "${NUM_JOBS}" \
         --joblog logs/parallel_consensus.log \
         -a "${cmd_file}" || {
           echo -e "\n[!] Warning: Some jobs in Step 3 failed."
           echo "Check 'logs/parallel_consensus.log' for details."
         }

# Count failures for a quick summary
FAILED_COUNT=$(awk '$7 != 0 {count++} END {print count+0}' logs/parallel_consensus.log)

# some of these jobs may fail if the consensus sequences are too short or there are not enough of them
if [ "$FAILED_COUNT" -gt 0 ]; then
    echo "Summary: $FAILED_COUNT task(s) failed. Proceeding to Step 4 with successful outputs..."
else
    echo "Parallel processing finished successfully with no errors."
fi
# ───────────────────────────────────────────────
# Step 4: Median K-mer Counting (R)
# ───────────────────────────────────────────────
echo -e "\nRunning Step 4: Median K-mer Counting..."

# remove previous file to avoid further concatenation
rm -f consensus/highest_median_kmer_samples.csv

Rscript scripts/median_kmer_counting.R "$@"

# ───────────────────────────────────────────────
# Step 5: Multi-Median Length Plotting (R)
# ───────────────────────────────────────────────
echo -e "\nRunning Step 5: Generating Median Length Plots..."
Rscript scripts/multi_median_length.R

# ───────────────────────────────────────────────
# Step 6: Concatenate Files (R)
# ───────────────────────────────────────────────
echo -e "\nRunning Step 6: Concatenating FASTA files..."
Rscript scripts/concatenate_files.R

echo -e "\nAll processing steps complete!"