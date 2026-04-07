#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo "Error: You must specify at least one experiment (e.g. oligo_1 oligo_2 oligo_3 streptavidin)"
  echo "Usage: $0 <exp1> <exp2> ..."
  exit 1
fi

mkdir -p cmds logs bait_lib_kmers data/fastas

echo "Generating labeling commands..."
: > cmds/label_kmers.sh

for exp in "$@"; do
  mkdir -p "results/${exp}/labeled"

  for k in {5..15}; do
    shopt -s nullglob
    for file in results/${exp}/enrichment_nolib/*_${k}mers_*.csv; do

      # pick bait file: exp-specific if present, otherwise universal
      bait_file="bait_lib_kmers/${exp}_comp_${k}mers.csv"
      if [ ! -f "$bait_file" ]; then
        bait_file="bait_lib_kmers/universal_comp_${k}mers.csv"
      fi

      # (optional) if even universal is missing, skip this one
      if [ ! -f "$bait_file" ]; then
        echo "Skipping: no bait file for ${exp}, k=${k} (missing exp-specific and universal)"
        continue
      fi

      fname="$(basename "${file}" .csv)"

      echo "python3 ./scripts/label_kmers.py \
        --results_file ${file} \
        --lib_file bait_lib_kmers/library_${k}mers.csv \
        --bait_file ${bait_file} \
        --outfile results/${exp}/labeled/${fname}_labeled.csv" \
      >> cmds/label_kmers.sh

    done
    shopt -u nullglob
  done
done


echo "Running labeling commands..."
bash cmds/label_kmers.sh
