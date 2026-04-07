# exit on any error
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

# make sure the appropriate directories are present
mkdir -p cmds logs results/${name_of_substrate}/counts

# memory intensive, takes awhile (hours to days: generating a separate cmds file for re-doing only streptavidin)
for k in {5..15}; do
    for infile in `ls results/${name_of_substrate}/graphs/*${k}mers.csv`; do
    name="$(basename $infile .csv)"
    outfile="${name}_counts.csv"
    echo "python3 ./scripts/collapse_kmer_uniques.py \
    --input_file ${infile} \
    --outfile results/${name_of_substrate}/counts/${outfile} \
    | tee -a logs/${name}_counts.log"
    done
done > cmds/${name_of_substrate}_generate_kmer_counters.sh

echo "Generating k-mer counting commands for ${name_of_substrate}..."

cat cmds/${name_of_substrate}_generate_kmer_counters.sh | parallel --bar --eta -j4

echo "Kmers counted successfully ♡"