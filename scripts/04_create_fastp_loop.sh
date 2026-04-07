#!/usr/bin/env bash

# Exit immediately on error
set -e

# create the necessary directories
mkdir -p data/qc/fastp

# define variables
name_of_substrate=$1

while IFS=$'\t' read -r -a y; do #read the R1R2 file
    o1="$(basename ${y[0]} .fastq.gz)"
	o2="$(basename ${y[1]} .fastq.gz)"
    merg="${o1%_L00*}.merged.fpqc.fastq.gz" # this is a new line
	echo "fastp -i ${y[0]} -I ${y[1]} \
	-o data/qc/${name_of_substrate}/${o1}.fpqc.fastq.gz \
	-O data/qc/${name_of_substrate}/${o2}.fpqc.fastq.gz \
	--failed_out data/qc/fastp/failed/${o1}.failed.txt \
	--merge --merged_out data/merged/${name_of_substrate}/${merg} \
	--length_limit 100 \
	--json data/qc/fastp/${name_of_substrate}_${o1}.json \
	--html data/qc/fastp/${name_of_substrate}_${o1}.html"
done < data/meta/raw_read_paths/${name_of_substrate}_R1R2.txt > cmds/fastp_${name_of_substrate}.cmds
echo "Fastp command loop for ${name_of_substrate} was created successfully ♡ "