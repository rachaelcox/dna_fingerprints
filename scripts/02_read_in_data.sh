#!/usr/bin/env bash

# check that two arguments are given
if [ $# -lt 2 ]; then
  echo "Usage: $0 <name_of_substrate> <source_directory>"
  exit 1
fi

# get variables
name_of_substrate=$1
source_dir=$2

# make sure the source directory exists
if [ ! -d "$source_dir" ]; then
  echo "Source directory not found: $source_dir"
  exit 1
fi

# create a variable to store the matched files
matched_files=()

# if the files exist then put them into matched_files; include symbolic links
for file in "$source_dir"/*; do
  # skip if no matches
  [ -e "$file" ] || [ -L "$file" ] || continue
  matched_files+=("$file")
done

# define the target directory relative to current project and make directory if it doent exist
target_dir="data/raw/${name_of_substrate}"
mkdir -p "$target_dir"

# if files are not found then throw an error
if [ ${#matched_files[@]} -eq 0 ]; then
  echo "No files or symlinks found in directory: $source_dir"
  exit 1
fi

# for the existing files
echo "${#matched_files[@]} files found in $source_dir"
echo "Linking ${#matched_files[@]} files to $target_dir ..."


for file in "${matched_files[@]}"; do
    if [ -L "$file" ] || [ -e "$file" ]; then
        ln -s "$file" "$target_dir/"
        echo "   → Linked $(basename "$file")"
    else
        echo "   Skipping missing file: $file"
    fi
done

echo "Linking complete for ${#matched_files[@]} files for ${name_of_substrate} ♡"