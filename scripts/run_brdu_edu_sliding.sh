#!/bin/bash
#SBATCH --job-name=run_brdu_edu_sliding
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=serial

# This shell script is to run 'brdu_edu_sliding_windows_from_bam.py' on an HPC cluster

# === Loading modules ===
module purge
module load slurm python pandas matplotlib

# === Paths and parameters ===
FILTERED_BED="/s/sansam-lab/path/to/filtered_bedfile"
BRDU_PREFIX="/s/sansam-lab/path/to/directory/to_store/All_bedgraphs/BrdU_"
EDU_PREFIX="/s/sansam-lab/path/to/directory/to_store/All_bedgraphs/EdU_"   #Keep the directories common for both the EdU and BrdU bedgraphs
BAM_FILE="/s/sansam-lab/path/to/directory/containing/the/corresponding_bamfile/"
WINDOW_SIZE=100
STEP_SIZE=10

# === Run the Python script ===
python /s/sansam-lab/path/to/Python_script/scripts/brdu_edu_sliding_windows_from_bam.py \
  --bam "$BAM_FILE" \
  --brdu_output "$BRDU_PREFIX" \
  --edu_output "$EDU_PREFIX" \
  --window_bp "$WINDOW_SIZE" \
  --step_size_bp "$STEP_SIZE" \
  --filtered "$FILTERED_BED" \
  --processes 16 \
  --output_per_read
  #--num_reads 2 \ 
  #--write_wig

