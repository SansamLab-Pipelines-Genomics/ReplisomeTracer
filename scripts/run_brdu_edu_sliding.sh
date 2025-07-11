#!/bin/bash
#SBATCH --job-name=run_brdu_edu_sliding
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=serial

# This shell script is to run 'brdu_edu_sliding_windows_from_bam.py' on a HPC cluster

# === Loading modules ===
module purge
module load slurm python pandas numpy matplotlib

# === Paths and parameters ===
FILTERED_BED="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/forkSense/Aux_1/final_filtered_trial.bed"
BRDU_PREFIX="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/analyze_bams/Aux_trial/BrdU_"
EDU_PREFIX="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/analyze_bams/Aux_trial/EdU_"
BAM_FILE="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/detect_combined/Aux_1.bam"
WINDOW_SIZE=100
STEP_SIZE=10

# === Run the script ===
python /s/sansam-lab/Try2_RepNano/RepNano/scripts/brdu_edu_sliding_windows_from_bam.py \
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

