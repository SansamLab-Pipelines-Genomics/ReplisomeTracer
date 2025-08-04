'''
Script Name: plot_summary.py

Description: This script processes the paired Edu-Brdu bedgraph files produced for all the 
reads (by the brdu_edu_sliding_windows_from_bam.py script output) to generate a smoothed 
aggregate plot of the Brdu vs Edu track distribution, and the corresponding aggregated dataframe
as a <.csv> file.

Usage: Run this python script using the shell script 'run_plot_summary.sh' The arguments 
to be passed will be in the shell script. They are amendable. Instructions on how to run 
it is in more detail within the header of the shell script.

Author: Krishanu Dhar

Date: July 10, 2025 (Last edited)  
'''

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import sem
from scipy.ndimage import uniform_filter1d

def read_bedgraph_tracks(directory):
    '''
    Reads all .bedgraph files in a directory and pairs BrdU and EdU files for the same sample.
    Steps:
    1. List all .bedgraph files in the directory.
    2. Load each file into a DataFrame.
    3. Extract unique sample basenames.
    4. Match BrdU__<sample> with EdU__<sample> and build a merged DataFrame with:
        (a) Chromosome, start, end
        (b) BrdU signal, EdU signal
        (c) diff = BrdU - EdU
    5. Returns: List of paired DataFrames.
    '''
    files = glob.glob(os.path.join(directory, "*.bedgraph"))
    dataframes = {os.path.basename(f): pd.read_csv(f, sep='\t', header=None) for f in files}
    basenames = [os.path.basename(f).split('__')[1] for f in files if '__' in os.path.basename(f)]
    unique_basenames = list(set(basenames))

    paired_tracks = []
    for base in unique_basenames:
        brdu_key = f"BrdU__{base}"
        edu_key = f"EdU__{base}"
        if brdu_key in dataframes and edu_key in dataframes:
            brdu = dataframes[brdu_key]
            edu = dataframes[edu_key]
            df = pd.DataFrame({
                "chromosome": brdu.iloc[:, 0],
                "start": brdu.iloc[:, 1],
                "end": brdu.iloc[:, 2],
                "BrdU": brdu.iloc[:, 3],
                "EdU": edu.iloc[:, 3],
            })
            df["diff"] = df["BrdU"] - df["EdU"]
            paired_tracks.append(df)
            #paired_tracks.append((base.replace(".bedgraph", ""), df))  # Use when running run_single_pair_boundary.sh
    return paired_tracks

def align_diff_by_minimum(df, window_size=500):
    '''
    Aligns tracks by centering the minimum of the smoothed difference signal.
    Steps:
    1. Smooth diff using a rolling uniform filter.
    2. Find where this smoothed signal hits min and max.
    3. If min appears after max, the data is reversed.
    4. Realign x-axis so that minimum is set to 0.
    5. Returns: The aligned DataFrame with x as the new relative position.
    '''
    df["diff_roll"] = pd.Series(uniform_filter1d(df["diff"].fillna(0), size=window_size, mode='nearest'))
    min_start = df.loc[df["diff_roll"].idxmin(), "start"]
    max_start = df.loc[df["diff_roll"].idxmax(), "start"]
    if min_start > max_start:
        df["start"] = df["start"].iloc[::-1].reset_index(drop=True)
    min_start = df.loc[df["diff_roll"].idxmin(), "start"]
    df["x"] = df["start"] - min_start
    return df

def summarize_combined_tracks(combined_df):
    '''
    Aggregates the aligned data from all samples.
    Steps:
    1. Groups rows by x (relative position).
    2. Computes:
        median_diff, mean, std, count
    3. Adds:
        se: standard error = std / sqrt(n)
        ci_lower and ci_upper: 95% confidence interval
    4. Returns: Summarized DataFrame.
    '''
    summary = combined_df.groupby("x")["diff"].agg(
        median_diff="median",
        mean="mean",
        sd="std",
        n="count"
    ).reset_index()
    summary["se"] = summary["sd"] / np.sqrt(summary["n"])
    summary["ci_lower"] = summary["mean"] - 1.96 * summary["se"]
    summary["ci_upper"] = summary["mean"] + 1.96 * summary["se"]
    return summary

def detect_signal_boundaries(summary_df):
    '''
    Detects domain boundaries based on median_diff.
    Steps:
    1. Smooth median_diff again.
    2. Find first location where signal drops below -0.1 → edu_start
    3. Then find first point > 0 after edu_start → brdu_start
    4. Then, 100 units later, find where it drops below 0.01 → brdu_end
    5. Returns: Same DataFrame with boundary columns added.
    '''
    roll_median = pd.Series(uniform_filter1d(summary_df["median_diff"].fillna(0), size=101, mode='nearest'))
    edu_start_idx = roll_median[roll_median < -0.1].index
    summary_df = summary_df.reset_index(drop=True)

    edu_start = summary_df.iloc[edu_start_idx[0]]["x"] if not edu_start_idx.empty else np.nan

    brdu_start = np.nan
    if not np.isnan(edu_start):
        brdu_candidates = summary_df[summary_df["x"] >= edu_start]
        brdu_above = brdu_candidates[brdu_candidates["median_diff"] > 0]
        if not brdu_above.empty:
            brdu_start = brdu_above.iloc[0]["x"]

    brdu_end = np.nan
    if not np.isnan(brdu_start):
        later_points = summary_df[summary_df["x"] > brdu_start + 100]
        roll_median_later = pd.Series(uniform_filter1d(later_points["median_diff"].fillna(0), size=101, mode='nearest'))
        brdu_end_candidates = roll_median_later[roll_median_later <= 0.01].index
        if not brdu_end_candidates.empty:
            brdu_end = later_points.iloc[brdu_end_candidates[0]]["x"]

    summary_df["edu_start"] = edu_start
    summary_df["brdu_start"] = brdu_start
    summary_df["brdu_end"] = brdu_end
    return summary_df

def make_aggregate_plot_df(bedgraph_directory, sample_name):
    '''
    High-level function to process all steps:
    Steps:
    1. Read and pair tracks.
    2. Align each track.
    3. Add sample ID to each.
    4. Concatenate them.
    5. Summarize into one DataFrame.
    6. Clip to -40kb to +60kb window.
    7. Add sample name and detect boundaries.
    8. Returns: Final summary DataFrame.
    '''
    replication_tracks = read_bedgraph_tracks(bedgraph_directory)
    aligned = [align_diff_by_minimum(df) for df in replication_tracks]
    for i, df in enumerate(aligned):
        df["id"] = i
    combined_df = pd.concat(aligned, ignore_index=True)
    summary_df = summarize_combined_tracks(combined_df)
    summary_df = summary_df[(summary_df["x"] >= -40000) & (summary_df["x"] <= 60000)]
    summary_df["sample_name"] = sample_name
    summary_df = detect_signal_boundaries(summary_df)
    return summary_df

def plot_domain_boundaries(summary_df, output_prefix):
    ''' 
    Extract boundary positions (assumes they're constant across the df), and 
    plots the median_diff over x, marking:
    edu_start, brdu_start and brdu_end
    Saves as .png, .pdf, and .svg.
    '''
    edu = summary_df['edu_start'].iloc[0]
    brdu = summary_df['brdu_start'].iloc[0]
    end = summary_df['brdu_end'].iloc[0]

    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(summary_df['x'], summary_df['median_diff'], color='red', linewidth=1)

    # Draw and label each boundary
    for pos, name in [(edu, 'edu_start'),
                      (brdu, 'brdu_start'),
                      (end, 'brdu_end')]:
        ax.axvline(pos, linestyle='--', linewidth=1)
        # place the text just above the bottom of the plot:
        ymin, ymax = ax.get_ylim()
        ax.text(pos, ymin + 0.05*(ymax-ymin),
                f"{int(pos)}",
                rotation=0,
                va='bottom', ha='center',
                fontsize=8)

    ax.set_xlim(summary_df['x'].min(), summary_df['x'].max())
    ax.set_xlabel('Relative Position (x)')
    ax.set_ylabel('Median diff')
    ax.set_title(f"Signal with Domain Boundaries\n{output_prefix}")

    for ext in ['png','pdf','svg']:
        fig.savefig(f"{output_prefix}.{ext}", dpi=300, bbox_inches='tight')
    plt.close(fig)


