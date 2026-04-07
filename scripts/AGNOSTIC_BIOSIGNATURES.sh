#!/bin/bash

# Ensure the script stops if any command fails
set -e

echo ">>> Starting Agnostic biosignature discovery pipeline"

# 1. READ IN THE FILES AND SET UP DIRECTORIES
echo ">>> Step 1: Setting up directories..."
./scripts/space_aptamers_01.sh oligo1 /stor/work/Ellington/sanchita/oligos/data/raw/oligo_1_redux
./scripts/space_aptamers_01.sh oligo2 /stor/work/Ellington/sanchita/oligos/data/raw/oligo_2_redux
./scripts/space_aptamers_01.sh oligo3 /stor/work/Ellington/sanchita/oligos/data/raw/oligo_3_redux
./scripts/space_aptamers_01.sh oligo4 /stor/work/Ellington/sanchita/oligos/data/raw/oligo_b_redux
./scripts/space_aptamers_01.sh streptavidin /stor/work/Ellington/sanchita/oligos/data/raw/streptavidin_redux
./scripts/space_aptamers_01.sh Chitin_Only /stor/work/Ellington/sanchita/chitin_PFAO/data/raw/Chitin_Only
./scripts/space_aptamers_01.sh glass161 /stor/work/Ellington/sanchita/glass161/data/raw/glass161

# 2. QC THE DATA
echo ">>> Step 2: Quality Control..."
./scripts/space_aptamers_02.sh oligo1 bg 0w 1w 3w
./scripts/space_aptamers_02.sh oligo2 bg 0w 1w 3w
./scripts/space_aptamers_02.sh oligo3 bg 0w 1w 3w
./scripts/space_aptamers_02.sh oligo4 bg 0w 1w 3w
./scripts/space_aptamers_02.sh streptavidin bg 0w 1w 3w
./scripts/space_aptamers_02.sh Chitin_Only bg 0w 1w 3w
./scripts/space_aptamers_02.sh glass161 bg 0w 1w 3w

# 3. GENERATE KMERS FOR SUBSTRATES
echo ">>> Step 3a: Generating Kmers..."
./scripts/space_aptamers_03a.sh degenerate_library oligo1 oligo2 oligo3 oligo4 streptavidin Chitin_Only glass161

# 4. MAKE DEBRUIJN GRAPHS
echo ">>> Step 3: Making De Bruijn Graphs..."
export NUM_JOBS=12
./scripts/space_aptamers_03.sh oligo1
./scripts/space_aptamers_03.sh oligo2
./scripts/space_aptamers_03.sh oligo3
./scripts/space_aptamers_03.sh oligo4
./scripts/space_aptamers_03.sh streptavidin
./scripts/space_aptamers_03.sh glass161
./scripts/space_aptamers_03.sh Chitin_Only

# 5. RUN DEGUST
echo ">>> Step 4: Running DEGUST..."
export NUM_JOBS=12
./scripts/space_aptamers_04.sh --conditions 0w,1w,3w oligo1 oligo2 oligo3 oligo4 streptavidin Chitin_Only glass161

# 6. LABEL THE ENRICHED KMER LIST
echo ">>> Step 4b: Labeling Enriched Kmers..."
./scripts/space_aptamers_04b.sh oligo1 oligo2 oligo3 oligo4 streptavidin Chitin_Only glass161

# 7. ANALYZE KMER SPACE
echo ">>> Step 4c: Analyze kmer space..."
export NUM_JOBS=12
./scripts/space_aptamers_04c.sh

# 8. FIGURE GENERATION
echo ">>> Step 5: Generating Figures..."
./scripts/space_aptamers_05.sh oligo1 oligo2 oligo3 oligo4 streptavidin Chitin_Only glass161

# 9. GENERATE THE CONSENSUS SEQUENCE
echo ">>> Step 6: Generating Consensus Sequences..."
./scripts/space_aptamers_06.sh oligo1 oligo2 oligo3 oligo4 streptavidin Chitin_Only glass161

# 10. BLAST ANALYSIS
echo ">>> Step 7: Running BLAST..."
./scripts/optional_space_aptamers_07.sh oligo1 oligo2 oligo3 oligo4

# 10. CLUSTERING ANALYSIS
echo ">>> Step 8: Final Clustering Analysis..."
./scripts/optional_space_aptamers_08.sh

echo ">>> Pipeline Complete!"