#!/usr/bin/env bash

# Define the substrate name (passed as first argument)
name_of_substrate=$1

# Make a folder to store the path lists
mkdir -p data/meta/raw_read_paths/

# Grab R1 and R2 file paths from the raw data directory
ls -1 data/raw/${name_of_substrate}/*R1*.fastq.gz > data/meta/raw_read_paths/${name_of_substrate}_R1.txt
ls -1 data/raw/${name_of_substrate}/*R2*.fastq.gz > data/meta/raw_read_paths/${name_of_substrate}_R2.txt

# make sure that each read as a pair
wc -l data/meta/raw_read_paths/*txt

# combine them into a single text file
paste data/meta/raw_read_paths/${name_of_substrate}_R1.txt data/meta/raw_read_paths/${name_of_substrate}_R2.txt > data/meta/raw_read_paths/${name_of_substrate}_R1R2.txt