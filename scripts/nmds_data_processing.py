#!/usr/bin/env python3
import pandas as pd
import os
import sys

def main():
    experiment_names = sys.argv[1:]
    if not experiment_names:
        print("Usage: ./nmds_process.py name1 name2 name3 ...")
        return

    input_dir = "figures/combined"
    nmds_data = []

    print(f"--- Step 1: Loading & Deduplicating ---")
    for name in experiment_names:
        file_path = os.path.join(input_dir, f"{name}_data.csv") # read in the data files after the combination
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            
            # grab the necessary columns
            subset = df[['edge', '0w', '1w', '3w']].copy() 
            
            # deduplicate the subsets
            subset = subset.drop_duplicates(subset='edge')
            
            # Rename columns
            subset.columns = ['edge', f"{name}_0w", f"{name}_1w", f"{name}_3w"]
            
            # append the data frames into a list
            nmds_data.append(subset)
            print(f"  [+] {name}: {subset.shape}")
        else:
            print(f"  [!] Skipping: {file_path} (not found)")

    if not nmds_data: # if there is no nmds data then quit
        print("No data loaded. Exiting.")
        return

    # outer merge of the nmds data
    print(f"\n--- Step 2: Outer Merging ---")
    nmds_full_merge = nmds_data[0]
    for i in range(1, len(nmds_data)): # merge the data frames on the edge and go through the list of the data frames
        nmds_full_merge = pd.merge(nmds_full_merge, nmds_data[i], on='edge', how='outer') 
        print(f"  Merged file {i+1}, current rows: {len(nmds_full_merge)}")

    # Save the raw merge just in case
    nmds_full_merge.to_csv("figures/combined/NMDS_FULL.csv", index=False)

    # Get rid of the negative data
    print(f"\n--- Step 3: Masking & Filtering ---")
    nmds_full_positive_only = nmds_full_merge.copy()
    
    # Get only the numeric columns for masking/dropping
    numeric_cols = nmds_full_positive_only.select_dtypes(include=['float64', 'int64']).columns
    
    # Replace negatives with NaN
    nmds_full_positive_only[numeric_cols] = nmds_full_positive_only[numeric_cols].mask(nmds_full_positive_only[numeric_cols] < 0)
    
    # Drop rows that are now ALL NaN in the numeric columns
    nmds_full_positive_only = nmds_full_positive_only.dropna(how='all', subset=numeric_cols)
    print(f"  Final Cleaned Shape: {nmds_full_positive_only.shape}")

    # slice data frame by their washes, this is hardcoded into 3, but can be changed if needed
    dff = nmds_full_positive_only
    dff_0 = dff[['edge'] + [c for c in dff.columns if c.endswith('_0w')]]
    dff_1 = dff[['edge'] + [c for c in dff.columns if c.endswith('_1w')]]
    dff_3 = dff[['edge'] + [c for c in dff.columns if c.endswith('_3w')]]

    print(f"  [Shapes] Full: {dff.shape}, 0w: {dff_0.shape}, 1w: {dff_1.shape}, 3w: {dff_3.shape}")



if __name__ == "__main__":
    main()