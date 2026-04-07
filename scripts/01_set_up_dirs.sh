
#!/usr/bin/env bash

# ensure the arguments are being input correctly
if [ $# -lt 1 ]; then
  echo "Usage: $0 <name_of_substrate>"
  exit 1
fi

# set variable for your substrate or dataset

name_of_substrate=$1

# make the directories needed
mkdir -p {data,cmds,logs,results,figures,scripts,bait_lib_kmers,consensus,blast}
mkdir -p "data/raw/${name_of_substrate}"
mkdir -p "data/fastas"
mkdir -p "data/kmer_space"
mkdir -p "data/qc/${name_of_substrate}"
mkdir -p "data/merged/${name_of_substrate}"
mkdir -p "data/processed/${name_of_substrate}"
mkdir -p "results/${name_of_substrate}"/{graphs,counts,enrichment,labeled}
mkdir -p figures/{combined,clean_tables,images}
mkdir -p consensus/{median_lengths,clustering}


# show completion
echo "Directory structure created for ${name_of_substrate}"