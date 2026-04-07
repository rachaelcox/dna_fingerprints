#!/usr/bin/env bash
set -e  # Exit on any error

# ───────────────────────────────
# Step 0:  Input checks
# ───────────────────────────────
if [ $# -lt 2 ]; then
  echo "Usage: $0 <name_of_substrate> <source_dir>"
  exit 1
fi

name_of_substrate=$1
source_dir=$2

echo "Running pipeline for: ${name_of_substrate}"
echo "Source directory: ${source_dir}"

if [ ! -d "$source_dir" ]; then
  echo "Source directory not found: $source_dir"
  exit 1
fi

# ───────────────────────────────
# Step 1: Set up Directories
# ───────────────────────────────
echo "Setting up directories..."
bash ./scripts/01_set_up_dirs.sh "${name_of_substrate}" || {
  echo "Step 1: Directory set up failed"
  exit 1
}

# ───────────────────────────────
# Step 2: Create symbolic links
# ───────────────────────────────
echo "Setting up dirctories..."
bash ./scripts/02_read_in_data.sh "${name_of_substrate}" "${source_dir}" || {
  echo "Step 2: Symbolic links could not be created"
  exit 1
}

# ───────────────────────────────
# Step 3: Create symbolic links
# ───────────────────────────────
echo "Identifying file paths..."
bash ./scripts/03_grab_fp.sh "${name_of_substrate}" || {
  echo "Step 3: File paths could not be reconciled"
  exit 1
}

# ───────────────────────────────
# Step 4: Generate fastp commands
# ───────────────────────────────
echo "Generating fastp commands..."
bash ./scripts/04_create_fastp_loop.sh "${name_of_substrate}" || {
  echo "Step 4 failed — could not create command file"
  exit 1
}


# ───────────────────────────────
# Done
# ───────────────────────────────
echo "Pre-processing complete ♡ "
echo "Substrate: ${name_of_substrate}"
echo "Now ready for QC"
