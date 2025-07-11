"""
Utilities for processing BAM files with DNA modifications. The parent python script
for this utilities script is 'brdu_edu_sliding_windows_from_bam.py'. 

This module provides helper functions for parsing DNA modification tags, computing 
modification frequencies using a sliding window approach, converting relative 
positions to genomic coordinates, and exporting results in BEDGRAPH or WIG format.

Functions:
----------
- parse_ml_mm(MM, ML)
- extract_read_data(read, window_size, step_size)
- compute_sliding_windows(mod_positions, window_size=100, step_size=10)
- convert_relative_to_abs_positions(reference_name, reference_start, reference_end, is_reverse, window_results)
- write_bedgraphs(brdu_output_file, edu_output_file, new_windows)
- write_wig(wig_output_file, chrom, start, values)

Author:
-------
Chris Sansam, Edited by Krishanu Dhar
Date:
-----
February 2025, Edited last in May 2025
"""

def parse_ml_mm(MM, ML):
    mod_positions = {}
    ml_probs = iter(ML)
    mm_entries = [entry for entry in MM.split(";") if entry]

    for entry in mm_entries:
        parts = entry.split(",")
        mod_type = parts[0]
        rel_positions = list(map(int, parts[1:]))
        abs_position = 0

        for rel_pos in rel_positions:
            abs_position += rel_pos + 1
            prob = next(ml_probs) / 255.0

            if abs_position not in mod_positions:
                mod_positions[abs_position] = {"BrdU": 0, "EdU": 0, "None": 1}

            if mod_type == "N+b?":
                mod_positions[abs_position]["BrdU"] = prob
            elif mod_type == "N+e?":
                mod_positions[abs_position]["EdU"] = prob

    for pos in mod_positions:
        total_mod_prob = mod_positions[pos]["BrdU"] + mod_positions[pos]["EdU"]
        mod_positions[pos]["None"] = max(0, 1 - total_mod_prob)

    return mod_positions


def extract_read_data(read, window_size, step_size):
    if not read.has_tag("MM") or not read.has_tag("ML"):
        return None

    return (
        read.reference_name,
        read.reference_start,
        read.reference_end,
        read.is_reverse,
        read.get_tag("MM"),
        read.get_tag("ML"),
        window_size,
        step_size
    )


def compute_sliding_windows(mod_positions, window_size, step_size):
    sorted_positions = sorted(mod_positions.keys())
    if not sorted_positions:
        return []

    min_pos = sorted_positions[0]
    max_pos = sorted_positions[-1]
    window_results = []

    for window_start in range(min_pos, max_pos - window_size + 1, step_size):
        window_end = window_start + window_size
        brdu_sum = 0.0
        edu_sum = 0.0
        t_count = 0

        for pos in range(window_start, window_end):
            if pos in mod_positions:
                brdu_sum += mod_positions[pos].get("BrdU", 0)
                edu_sum += mod_positions[pos].get("EdU", 0)
                t_count += 1

        if t_count > 0:
            brdu_freq = brdu_sum / t_count
            edu_freq = edu_sum / t_count
        else:
            brdu_freq = 0.0
            edu_freq = 0.0

        window_results.append((window_start, window_end, brdu_freq, edu_freq))

    return window_results


def convert_relative_to_abs_positions(reference_name, reference_start, reference_end, is_reverse, window_results):
    pos_list = tuple(int(tup[0]) for tup in window_results)
    window_size = window_results[0][1] - window_results[0][0]

    if not is_reverse:
        starts = tuple(reference_start + 1 + x for x in pos_list)
    else:
        starts = tuple(reference_end - 1 - x for x in pos_list)

    new_windows = tuple(
        (reference_name, x, x + window_size, window_results[i][2], window_results[i][3])
        for i, x in enumerate(starts)
    )

    return new_windows


def write_bedgraphs(brdu_output_file, edu_output_file, new_windows):
    with open(brdu_output_file, "w") as f:
        for window in new_windows:
            chrom, start, end, brdu_freq, _ = window
            f.write(f"{chrom}\t{start}\t{end}\t{brdu_freq:.6f}\n")
    print(f"BEDGRAPH file saved: {brdu_output_file}")

    with open(edu_output_file, "w") as f:
        for window in new_windows:
            chrom, start, end, _, edu_freq = window
            f.write(f"{chrom}\t{start}\t{end}\t{edu_freq:.6f}\n")
    print(f"BEDGRAPH file saved: {edu_output_file}")


def write_wig(wig_output_file, chrom, start, values):
    with open(wig_output_file, "w") as f:
        f.write(f"fixedStep chrom={chrom} start={start} step=1\n")
        for val in values:
            f.write(f"{val:.6f}\n")
    print(f"WIG file saved: {wig_output_file}")
