'''
Script Name: filtered_bedfiles.py

Description: This python scripts takes the bedfile outputs from forkSense, merges them and
filters them to generate a new bedfile. The resulting bedfiles will be fed into the 
'brdu_edu_sliding_windows_from_bam' Python script next. This script will still preserve the 
original metadata associated with each read ID from the original bedfiles.

Usage: Use the shell script 'run_filtered_bedfiles.sh' to run this Python script.

Author: Krishanu Dhar

Date: July 11, 2025 (Last edited)

'''
import pysam
import pandas as pd
import argparse

def get_merged_forks_filtered(left_path, right_path, origin_path, termination_path, read_id_col=3, score_col=8):
    left_df = pd.read_csv(left_path, sep=r"\s+", header=None, comment="#")
    right_df = pd.read_csv(right_path, sep=r"\s+", header=None, comment="#")
    origin_df = pd.read_csv(origin_path, sep=r"\s+", header=None, comment="#")
    termination_df = pd.read_csv(termination_path, sep=r"\s+", header=None, comment="#")

    merged_df = pd.concat([left_df, right_df], ignore_index=True)
    merged_ids = set(merged_df[read_id_col].astype(str).str.strip())
    exclude_ids = set(pd.concat([origin_df, termination_df])[read_id_col].astype(str).str.strip())

    filtered_ids = merged_ids - exclude_ids
    filtered_df = merged_df[merged_df[read_id_col].isin(filtered_ids)]
    #filtered_df = filtered_df[filtered_df[score_col] != -3.000000] # For excluding read IDs which have forks trailing off the ends of a read

    return filtered_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process BAM file and generate BrdU & EdU BEDGRAPH files."
    )
    parser.add_argument("--left", required=True, help="Path to left fork bedfile")
    parser.add_argument("--right", required=True, help="Path to right fork bedfile")
    parser.add_argument("--origin", required=True, help="Path to origin bedfile")
    parser.add_argument("--termination", required=True, help="Path to termination bedfile")
    parser.add_argument("--bed_output", required=True, help="Path to filtered bed file")
    args = parser.parse_args()
    bed_df = get_merged_forks_filtered(
		args.left,
		args.right,
		args.origin,
		args.termination)
    bed_df.columns = ['chrom', 'start', 'end', 'read_id', 'col4', 'col5', 'strand', 'length', 'score']
    bed_df.to_csv(args.bed_output, sep=" ", header=False, index=False)
    print(f"Filtered BED file saved to: {args.bed_output}")



