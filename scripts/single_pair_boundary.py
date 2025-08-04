#!/usr/bin/env python3
'''
Name: single_pair_boundary.py
To run: Run it with 'sbatch run_single_pair_boundary.sh'
'''
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from plot_summary import read_bedgraph_tracks, align_diff_by_minimum


def plot_each_track(input_dir):
    paired_tracks = read_bedgraph_tracks(input_dir)
    output_dir = os.path.join(input_dir, "annotated_track_plots")
    os.makedirs(output_dir, exist_ok=True)

    results = []
    check = 0

    for i, (sample_id, df) in enumerate(paired_tracks):
        #check += 1
        #if check == 21: break
        # if sample_id != "3bebeea1-f150-4853-a311-1f295dfd35c8":
        #     continue
        # before_aligned_path = os.path.join(output_dir, "before_3bebeea1-f150-4853-a311-1f295dfd35c8.bed")
        # df.to_csv(before_aligned_path, sep=" ", index=False)
        aligned_df = align_diff_by_minimum(df)
        # aligned_path = os.path.join(output_dir, "error_3bebeea1-f150-4853-a311-1f295dfd35c8.bed")
        # aligned_df.to_csv(aligned_path, sep=" ", index=False)
        # print(f"[✓] Exported: {aligned_path}")
        read_start = aligned_df["x"].min()
        read_end = aligned_df["x"].max()
        aligned_df = edu_boundaries_detect(aligned_df)

        edu_start = aligned_df['edu_start'].iloc[0]
        brdu_start = aligned_df['brdu_start'].iloc[0]
        brdu_end = aligned_df['brdu_end'].iloc[0]

        # --- PLOTTING ---
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(aligned_df["x"], aligned_df["diff"], label=sample_id, color='blue')
        ax.set_title(f"Aligned BrdU-EdU Track: {sample_id}")
        ax.set_xlabel("Relative Position (x)")
        ax.set_ylabel("BrdU - EdU (diff)")
        ax.axvline(0, linestyle="--", color="gray", linewidth=0.8)
        ax.legend()

        for pos, name, color in [(edu_start, 'edu_start', 'red'), (brdu_start, 'brdu_start', 'darkgreen'), (brdu_end, 'brdu_end', 'black')]:
            if not np.isnan(pos):
                ax.axvline(pos, linestyle='--', linewidth=1.5, color=color)
                ymin, ymax = ax.get_ylim()
                ax.text(pos, ymin + 0.05*(ymax - ymin), f"{int(pos)}",
                        rotation=0, va='bottom', ha='center', fontsize=8, color=color)

        plot_path = os.path.join(output_dir, f"{sample_id}1_plot.png")
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        # --- CSV ENTRY ---
        if (np.isnan(edu_start) or (edu_start - read_start) < 1500) or (np.isnan(brdu_end) or abs(brdu_end - read_end) < 1500):
            continue  # skip writing to CSV, but keep the plot

        edu_track = brdu_start - edu_start if not np.isnan(edu_start) and not np.isnan(brdu_start) else np.nan
        brdu_track = brdu_end - brdu_start if not np.isnan(brdu_end) and not np.isnan(brdu_start) else np.nan
        results.append({
            "read_id": sample_id,
            "edu_start": edu_start,
            "brdu_start": brdu_start,
            "brdu_end":brdu_end,
            "edu_track": edu_track,
            "brdu_track": brdu_track,
            "edu_speed": edu_track/10000,
            "brdu_speed": brdu_track/10000,
            "Fork_track": edu_track + brdu_track,
            "Fork_speed": (edu_track + brdu_track)/20000
        })

        # if sample_id == "3bebeea1-f150-4853-a311-1f295dfd35c8":
        #     break

    # --- WRITE CSV TO OUTPUT DIR ---
    results_df = pd.DataFrame(results)
    csv_path = os.path.join(output_dir, "edu_brdu_summary_new.csv")
    results_df.to_csv(csv_path, index=False)
    print(f"[✓] Exported: {csv_path}")


def edu_boundaries_detect(summary_df):
    # Step 1: Sort by x
    summary_df = summary_df.sort_values("x").reset_index(drop=True)

    smoothed_diff = pd.Series(
        uniform_filter1d(summary_df["diff"].fillna(0), size=101, mode='nearest')
    )
    summary_df["smoothed_diff"] = smoothed_diff

    # Step 3: Find EdU start
    edu_mask = (
        (summary_df["smoothed_diff"] < -0.1) &
        (summary_df["x"] < 0) &
        (summary_df["x"] > -10000)
    )
    edu_idx = edu_mask[edu_mask].index
    edu_start = summary_df.loc[edu_idx[0], "x"] if not edu_idx.empty else np.nan

    # Step 4: Find BrdU start
    brdu_start = np.nan
    if not np.isnan(edu_start):
        brdu_candidates = summary_df[(summary_df["x"] >= edu_start) & (summary_df["x"] > 0)]
        brdu_above = brdu_candidates[brdu_candidates["diff"] > 0]
        if not brdu_above.empty:
            brdu_start = brdu_above.iloc[0]["x"]

    # Step 5: BrdU end detection using smoothed signals
    brdu_end = np.nan
    if not np.isnan(brdu_start):
        df_filtered = summary_df.dropna(subset=["diff", "BrdU"]).copy()

        # Step 5.1: Smooth signals
        df_filtered["BrdU_smooth"] = df_filtered["BrdU"].rolling(window=201, center=True, min_periods=1).mean()
        df_filtered["diff_smooth"] = df_filtered["diff"].rolling(window=201, center=True, min_periods=1).mean()
        df_filtered["diff_div_brdu_smooth"] = df_filtered["diff_smooth"] / df_filtered["BrdU_smooth"]

        # Step 5.2: Window after brdu_start (15,000 bp)
        window_df = df_filtered[(df_filtered["x"] > brdu_start) & (df_filtered["x"] <= brdu_start + 15000)]

        if not window_df.empty:
            # Step 5.3: Find peak of smoothed BrdU
            brdu_peak_idx = window_df["BrdU_smooth"].idxmax()
            brdu_peak_x = df_filtered.loc[brdu_peak_idx, "x"]

            # Step 5.4: Search for drop zone from peak onward (on smoothed signals)
            drop_zone = df_filtered[df_filtered["x"] > brdu_peak_x]
            drop_zone = drop_zone[
                (drop_zone["diff_smooth"] <= 0.75) &
                (drop_zone["diff_div_brdu_smooth"] <= 0.65)
            ]

            if not drop_zone.empty:
                brdu_end = drop_zone.iloc[0]["x"]


    # Final output
    summary_df["edu_start"] = edu_start
    summary_df["brdu_start"] = brdu_start
    summary_df["brdu_end"] = brdu_end

    return summary_df


