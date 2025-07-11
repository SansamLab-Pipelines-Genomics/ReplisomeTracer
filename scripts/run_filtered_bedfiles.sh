#!/bin/bash
#SBATCH --job-name=run_filtered_bedfiles

# === Loading modules ===
module purge
module load slurm python pandas

# === Paths and parameters ===
LEFT_FORKS="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/forkSense/Aux_1/leftForks_DNAscent_forkSense.bed"
RIGHT_FORKS="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/forkSense/Aux_1/rightForks_DNAscent_forkSense.bed"
ORIGIN="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/forkSense/Aux_1/origins_DNAscent_forkSense.bed"
TERMINATION="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/forkSense/Aux_1/terminations_DNAscent_forkSense.bed"
BED_OUTPUT="/s/sansam-lab/CS_20250702_Rep/ReplisomeTracer/DNAscent_Snakemake/results/forkSense/Aux_1/final_filtered_trial.bed"

python /s/sansam-lab/Try2_RepNano/RepNano/scripts/filtered_bedfiles.py \
  --left "$LEFT_FORKS" \
  --right "$RIGHT_FORKS" \
  --origin "$ORIGIN" \
  --termination "$TERMINATION" \
  --bed_output "$BED_OUTPUT"