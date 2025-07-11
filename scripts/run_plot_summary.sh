#!/bin/bash
#SBATCH --cpus-per-task=16

# This shell script is to run the python script 'plot_summary.py'

# How to run:
# Give permission if not: chmod +x run_plot.sh
# Run this in the same directory the scripts are stored:
# Command: ./run_plot.sh <bedgraph_dir> <sample_name>
# The outputs will be placed in the same directory

# ===Loading modules===
module purge
ml slurm python pandas python matplotlib scipy

set -e

BEDGRAPH_DIR="$1"
SAMPLE_NAME="$2"


# Run Python analysis
python3 -c "
from plot_summary import make_aggregate_plot_df, plot_domain_boundaries
df = make_aggregate_plot_df('$BEDGRAPH_DIR', '$SAMPLE_NAME')
df.to_csv('${SAMPLE_NAME}_summary.csv', index=False)
plot_domain_boundaries(df, '${SAMPLE_NAME}_plot')
"

echo "Done. Outputs: ${SAMPLE_NAME}_summary.csv and ${SAMPLE_NAME}_plot.(png|pdf|svg)"
