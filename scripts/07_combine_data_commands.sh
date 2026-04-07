# exit on any error
set -e

# set a default cutoff
cutoff=20

# check if the first argument is the optional cutoff
if [ "$1" == "--cutoff" ]; then
  # ensure a value follows
  if [ -z "$2" ]; then
    echo "Error: --cutoff flag requires a numeric value"
    exit 1
  fi
  cutoff="$2"
  shift 2  # remove the first two arguments (--cutoff <value>)
fi

# check the arguments if cutoff is not there
if [ $# -lt 1 ]; then
  echo "Error: You must specify at least one substrate after the cutoff value"
  exit 1
fi

# ensure directories are present
mkdir -p cmds logs

# create a fresh command file at the start
cmd_file="cmds/combine_data_${cutoff}c.sh"
: > "$cmd_file"   # clears it once, before the loop

# read the arguments correctly
for dir in "$@"; do
  if [ ! -d "results/${dir}/counts" ]; then
    echo "Skipping ${dir} (no counts/ directory found)"
    continue
  fi

  mkdir -p "results/${dir}/counts_combined"

  for k in {5..15}; do
    echo "python3 ./scripts/combine_data.py \
      --data_dir results/${dir}/counts/ \
      --k ${k} \
      --cutoff ${cutoff} \
      --out_prefix results/${dir}/counts_combined/${dir}_${k}mers_${cutoff}c_combined" \
      >> "$cmd_file"
  done
done

echo "Running combine_data.py for substrates: $@ (cutoff=${cutoff})..."
cat "$cmd_file" | parallel --bar --eta -j4
echo "All combine_data jobs completed ♡"
