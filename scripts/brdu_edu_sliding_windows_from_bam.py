#!/usr/bin/env python3
"""
Script Name: brdu_edu_sliding_windows_from_bam.py

Description:
This script processes a BAM file to compute probability-weighted BrdU and EdU 
modification frequencies using a sliding window approach. This Python script takes the help of a 
utilities file 'slidingWindowsBedgraphs_utils.py' to process the read data and generate bedgraph
files for visualization. You can optionally write separate files per read and WIG files, too.
The filtered_bed file is simply a bed file containing the selected read IDs that one wants to process. 

Usage:
    Run this Python script with the shell script 'run_brdu_edu_sliding.sh'. The arguments provided in
    the shell script are 'default' arguments. Edit in your specifics.

Arguments:
    -b, --bam          Path to input BAM file (required) (Also, make sure the indexed BAM file is also in the same directory)
    -o1, --brdu_output Path to output BrdU BEDGRAPH file (or prefix if --output_per_read is used)
    -o2, --edu_output  Path to output EdU BEDGRAPH file (or prefix if --output_per_read is used)
    -p, --processes    Number of CPU cores to use (default: auto-detect)
    -n, --num_reads    Number of reads to process (default: all)
    -w, --window_bp    Window size in base pairs (default: 100)
    -s, --step_size_bp Step size in base pairs (default: 10)
    -r, --filtered     BED file with read IDs to filter
    --output_per_read  Output one BEDGRAPH, each for EdU and BrdU, per read
    --write_wig        Also write WIG files per read

Author:
    Chris Sansam, Edited by Krishanu Dhar

Date:
    February 2025, Edited Last on May 2025
"""

import pysam
import argparse
import pandas as pd
import sys
import os
import itertools
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
LIB_PATH = os.path.join(SCRIPT_DIR, "lib")
sys.path.insert(0, LIB_PATH)

import slidingWindowsBedgraphs_utils as sw


def process_read(read_data):
    if read_data is None:
        return []

    (
        reference_name,
        reference_start,
        reference_end,
        is_reverse,
        MM,
        ML,
        window_size,
        step_size,
        read_name
    ) = read_data

    mod_positions = sw.parse_ml_mm(MM, ML)
    window_results = sw.compute_sliding_windows(mod_positions, window_size, step_size)
    abs_results = sw.convert_relative_to_abs_positions(
        reference_name, reference_start, reference_end, is_reverse, window_results
    )
    return (read_name, abs_results)


def process_bam_and_write_bedgraphs(
    bam_path,
    filtered_bed,
    brdu_output,
    edu_output,
    window_size,
    step_size,
    num_processes=None,
    num_reads=None,
    output_per_read=False,
    write_wig=False
):
    if num_processes is None:
        num_processes = max(1, cpu_count() - 1)

    print(f"Using {num_processes} CPU cores for multiprocessing...")

    os.makedirs(os.path.dirname(brdu_output), exist_ok=True)
    os.makedirs(os.path.dirname(edu_output), exist_ok=True)

    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        if filtered_bed:
            bed_df = pd.read_csv(filtered_bed, sep=r"\s+", header=None)
            bed_df.columns = ['chrom', 'start', 'end', 'read_id', 'col4', 'col5', 'strand', 'length', 'score']
            read_ids_of_interest = set(bed_df['read_id'].astype(str).str.strip())
            read_chunks = [read for read in bamfile if read.query_name in read_ids_of_interest]
        else:
            read_chunks = list(bamfile.fetch())

        total_reads = len(read_chunks) if num_reads is None else min(len(read_chunks), num_reads)

        valid_reads = (
            read for read in read_chunks if sw.extract_read_data(read, window_size, step_size)
        )
        read_data_list = (
            tuple(list(sw.extract_read_data(read, window_size, step_size)) + [read.query_name])
            for read in itertools.islice(valid_reads, num_reads)
        )

        with Pool(processes=num_processes) as pool:
            results = list(
                tqdm(
                    pool.imap(process_read, read_data_list),
                    total=total_reads,
                    desc="Processing BAM Reads",
                    unit="read",
                )
            )

    if output_per_read:
        for read_name, windows in results:
            brdu_path = f"{brdu_output}_{read_name}.bedgraph"
            edu_path = f"{edu_output}_{read_name}.bedgraph"
            sw.write_bedgraphs(brdu_path, edu_path, windows)

            if write_wig:
                if windows:
                    chrom = windows[0][0]
                    start = windows[0][1]
                    brdu_vals = [w[3] for w in windows]
                    edu_vals = [w[4] for w in windows]
                    sw.write_wig(f"{brdu_output}_{read_name}.wig", chrom, start, brdu_vals)
                    sw.write_wig(f"{edu_output}_{read_name}.wig", chrom, start, edu_vals)
    else:
        all_window_results = [item for _, sublist in results for item in sublist]
        sw.write_bedgraphs(brdu_output, edu_output, all_window_results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process BAM file and generate BrdU & EdU BEDGRAPH files."
    )

    parser.add_argument("-b", "--bam", required=True, help="Path to input BAM file")
    parser.add_argument("-o1", "--brdu_output", required=True, help="Output BrdU BEDGRAPH file")
    parser.add_argument("-o2", "--edu_output", required=True, help="Output EdU BEDGRAPH file")
    parser.add_argument("-p", "--processes", type=int, default=None, help="Number of CPU cores")
    parser.add_argument("-n", "--num_reads", type=int, default=None, help="Number of reads to process")
    parser.add_argument("-w", "--window_bp", type=int, default=100, help="Window size in bp")
    parser.add_argument("-s", "--step_size_bp", type=int, default=10, help="Step size in bp")
    parser.add_argument("-r", "--filtered", default=None, help="Path to BED file with read IDs")
    parser.add_argument("--output_per_read", action="store_true", help="Output one BEDGRAPH per read")
    parser.add_argument("--write_wig", action="store_true", help="Also write WIG files per read")

    args = parser.parse_args()

    process_bam_and_write_bedgraphs(
        args.bam,
        args.filtered,
        args.brdu_output,
        args.edu_output,
        args.window_bp,
        args.step_size_bp,
        args.processes,
        args.num_reads,
        args.output_per_read,
        args.write_wig
    )
