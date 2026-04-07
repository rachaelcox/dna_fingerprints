# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# To run this the following files should exist for each substrate of interest (exp)
# complementary fastas: data/fastas/${exp}_comp.fasta
# direct match fastas: input_fasta data/fastas/${exp}.fasta

#!/usr/bin/env bash

# exit on any error
set -e

# ───────────────────────────────────────────────
# Step 0: Input check
# ───────────────────────────────────────────────
if [ $# -lt 2 ]; then
  echo "Usage: $0 <prey_name> <exp1> <exp2> ..."
  echo "Example: $0 universal oligo1 oligo2 oligo3"
  exit 1
fi

prey_name="$1"
shift

prey_fasta="data/fastas/${prey_name}_prey.fasta"
if [ ! -f "${prey_fasta}" ]; then
  echo "Error: prey FASTA not found: ${prey_fasta}"
  exit 1
fi


# ───────────────────────────────────────────────
# Step 1: Setup
# ───────────────────────────────────────────────
mkdir -p cmds logs bait_lib_kmers
mkdir -p data/fastas
mkdir -p bait_lib_kmers

# ───────────────────────────────────────────────
# step 2: generate complementary k-mers
# ───────────────────────────────────────────────
echo "generating complementary k-mers..."

universal_comp="data/fastas/universal_comp.fasta"

for k in {5..15}; do
  for exp in "$@"; do

    comp_fasta="data/fastas/${exp}_comp.fasta"

    if [ -f "${comp_fasta}" ]; then
      input_fasta="${comp_fasta}"
      echo "using ${comp_fasta}"
    elif [ -f "${universal_comp}" ]; then
      input_fasta="${universal_comp}"
      echo "${comp_fasta} not found, using universal_comp.fasta"
    else
      echo "error: no complementary fasta found for ${exp}"
      exit 1
    fi

    python3 ./scripts/generate_deBruijn_kmers.py \
      --input_fasta "${input_fasta}" \
      --kmer_size ${k} \
      --outfile bait_lib_kmers/${exp}_comp_${k}mers.csv

  done
done

echo "complementary k-mers generated"

# ───────────────────────────────────────────────
# Step 3: Generate library (prey) k-mers
# ───────────────────────────────────────────────
echo "Generating library (prey) k-mers..."
for k in {5..15}; do
  python3 ./scripts/generate_deBruijn_kmers.py \
    --input_fasta "${prey_fasta}" \
    --kmer_size ${k} \
    --outfile "bait_lib_kmers/library_${k}mers.csv"
done
echo "Library k-mers generated successfully"
