# Exit on any error
set -e

# set argument
name_of_substrate="$1"

# check that there are at least 2 arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <name_of_substrate>"
    echo "Example: $0 streptavidin"
    exit 1
fi

# default parallel job count (can override with NUM_JOBS)
num_jobs="${NUM_JOBS:-4}"

# make sure the correct directories exist
mkdir -p "results/${name_of_substrate}/graphs/" "cmds"

# check to see that there are FastQC/fastp outputs to analyze
input_pattern="data/qc/fastp/${name_of_substrate}*"
found_files=$(find data/qc/fastp -maxdepth 1 -type f -name "${name_of_substrate}*" | wc -l)

if [ "$found_files" -eq 0 ]; then
    echo " No matching files found for pattern: ${input_pattern}"
    echo " MultiQC cannot run without valid input files."
    exit 1
fi


# check to see if the infiles exist
if [ -z "$(ls data/processed/${name_of_substrate} 2>/dev/null)" ]; then
    echo "No input files found in data/processed/${name_of_substrate}"
    exit 1
fi


# low memory, quick generating a separate cmds file for re-doing only streptavidin

#run command loop
for k in {5..15}; do
    for infile in `ls data/processed/${name_of_substrate}`; do
    outfile=${infile%.merged*}_${k}mers.csv
    echo "python3 ./scripts/generate_deBruijn_kmers.py \
    --input_fasta data/processed/${name_of_substrate}/${infile} \
    --outfile results/${name_of_substrate}/graphs/${outfile} \
    --kmer_size ${k}"
    done
done > cmds/${name_of_substrate}_generate_dbgs.sh

#check progress
echo "Running $(wc -l < cmds/${name_of_substrate}_generate_dbgs.sh) graph generation jobs using $num_jobs cores..."

# run the actual command after generation
cat cmds/${name_of_substrate}_generate_dbgs.sh | parallel --bar -j${num_jobs}

# rename the files to make them more human readable
if ls ${name_of_substrate}/graphs/*.csv &>/dev/null; then
    for f in ${name_of_substrate}/graphs/*.csv; do
        newname=$(echo "$f" | sed 's/{\([0-9]\+\)}/\1/')
        mv "$f" "$newname"
    done
fi

echo "De Bruijn graphs successfully generated ♡"