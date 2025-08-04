#!/bin/bash

# ===Loading modules===
module purge
ml slurm python pandas matplotlib scipy
# Set variables
INPUT_DIR="/s/sansam-lab/20250731_Rep/ReplisomeTracer/analyze_bams/PARPi/"              # Change to your actual input folder

# Run the script
python3 -c "
from single_pair_boundary import plot_each_track, edu_boundaries_detect
plot_each_track('$INPUT_DIR')
"
echo "Done"
