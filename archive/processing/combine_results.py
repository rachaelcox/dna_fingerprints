"""
Script for combining k-mer enrichment results, optionally filtered at a given FDR threshold.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import time
import os
import re

''' Functions '''


''' Main '''
def main():

    # access arguments passed from the command line
    results_dir = args.results_dir
    exp_label = args.experiment_label
    fdr_threshold_str = args.fdr
    
    # convert FDR threshold to float if provided
    fdr_threshold = float(fdr_threshold_str) if fdr_threshold_str else None

    t0 = time.time()
    print('Reading in files ...')
    
    # get list of files
    all_files = os.listdir(results_dir)
    csv_files = [f for f in all_files if f.endswith('.csv')]
    
    # check if any CSV files were found
    if not csv_files:
        print(f"Error: No CSV files found in {results_dir}. Exiting.")
        return

    # read, process, and combine files
    list_of_dfs = []
    k_mer_regex = re.compile(r'(\d+)mer')
    
    for filename in csv_files:
        filepath = os.path.join(results_dir, filename)
        k_value = None
        match = k_mer_regex.search(filename)

        if match:
            # group 1 (the digits) is the k-mer size
            k_value = int(match.group(1))
            print(f"Detected k-mer size k={k_value} from filename: {filename}")
        else:
            print(f"Warning: Could not determine k-mer size from filename: {filename}. Skipping file.")
            continue # skip file if k cannot be determined
        
        # read the file into a df
        # the first column (k-mer sequence) is the index for a cleaner structure,
        try:
            current_df = pd.read_csv(filepath)
        except Exception as e:
            print(f"Warning: Could not read file {filepath}. Skipping. Error: {e}")
            continue
        
        # delete the 'confect' column (use .drop with axis=1 for column)
        if 'confect' in current_df.columns:
            current_df = current_df.drop(columns=['confect'])
        else:
            print(f"Warning: 'confect' column not found in {filename}. Skipping drop.")

        # add the new 'exp' column with the argument value
        current_df['exp'] = exp_label

        # add the new 'k' column with the extracted k-mer size
        current_df['k'] = k_value
        
        # optional: filter on the 'fdr' column
        if fdr_threshold is not None:
            # filter rows where the 'fdr' value is less than or equal to the threshold
            current_df = current_df[current_df['fdr'] <= fdr_threshold]
            print(f"Filtered {filename}: Kept {len(current_df)} rows below FDR {fdr_threshold*100}%.")
        
        list_of_dfs.append(current_df)
        print(f"Processed {filename} ({len(current_df)} rows).")

    # concatenate all dfs
    if list_of_dfs:
        df = pd.concat(list_of_dfs, ignore_index=True)
        print(f'Successfully combined {len(list_of_dfs)} files into one DataFrame with {len(df)} total rows.')
    else:
        print('No data to combine. Exiting.')
        return

    print('Writing out results ...')
    outfile = args.outfile
    # write the combined DataFrame to the specified output file
    df.to_csv(outfile, index=False)
    if args.write_edges:
        edge_outfile = outfile.replace('.csv', '_weighted_edges_only.csv')
        edge_df = df[['edge', '3w', 'k']]
        edge_df = edge_df.rename(columns={'3w':'weight'})
        edge_df.to_csv(edge_outfile, index=False)
    print('Done!')
    print(f'Total run time: {round((time.time()-t0)/60, 2)} minutes.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--results_dir", action="store", required=True,
                        help="(Required) Path to directory with enrichment results for multiple k.")

    parser.add_argument("--experiment_label", action="store", required=True,
                        help="(Required) Value to population the 'exp' column.")

    parser.add_argument("--outfile", action="store", required=True,
                        help="(Required) Outfile path/name (.csv)")

    parser.add_argument("--fdr", action="store", required=False,
                        help="(Optional) Filter for k-mers enriched below the specified FDR threshold, as a fraction (eg for 5% threshold, use 0.05)")

    parser.add_argument("--write_edges", action="store_true", required=False,
                        help="(Optional) Additionally, write out a file with edges and weights only, for each k.")

    # optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()