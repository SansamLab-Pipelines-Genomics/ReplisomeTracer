#!/bin/bash
#SBATCH --job-name=run_snakemake_pipeline
#Shell script to run the snakemake pipeline in your HPC cluster. You need to load the slurm module prior to running this script.

# Load required modules
module purge
module load slurm python pandas numpy matplotlib

# Run Snakemake with cluster support
snakemake -j 999 \
  --use-envmodules \
  --latency-wait 100 \
  --keep-incomplete \
  --cluster-config config/cluster_config.yml \
  --cluster "sbatch --account {cluster.account} \
                    --partition {cluster.partition} \
                    --cpus-per-task={cluster.cpus-per-task} \
                    --mem={cluster.mem} \
                    --output={cluster.output} \
                    --gres={cluster.gres} \
                    --constraint=\"{cluster.constraint}\""
