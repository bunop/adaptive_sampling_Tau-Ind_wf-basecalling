#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=96G
#SBATCH --error=tar_make.err
#SBATCH --output=tar_make.out
#SBATCH --job-name=tar_make
#SBATCH --account=IscrC_BioGPUPX    # account name
#SBATCH --partition=g100_usr_prod   # partition name (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)
#SBATCH --qos=normal                # quality of service (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)

# Check if conda or micromamba is available and activate the environment accordingly
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate R-4.5
elif command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook -s bash)"
    micromamba activate R-4.5
else
    echo "Neither conda nor micromamba found in PATH."
    exit 1
fi

# execute the tar_make command
if ! Rscript -e "targets::tar_make()"; then
    echo "Error: tar_make command failed."
    exit 1
fi

# prune the targets
if ! Rscript -e "targets::tar_prune()"; then
    echo "Error: tar_prune command failed."
    exit 1
fi
