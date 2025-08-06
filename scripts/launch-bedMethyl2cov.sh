#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00
#SBATCH --mem=16G
#SBATCH --error=bedMethyl2cov.err
#SBATCH --output=bedMethyl2cov.out
#SBATCH --job-name=bedMethyl2cov
#SBATCH --account=IscrC_BioGPUPX    # account name
#SBATCH --partition=g100_usr_prod   # partition name (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)
#SBATCH --qos=normal                # quality of service (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)

# Check if conda or micromamba is available and activate the environment accordingly
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate python3.12
elif command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook -s bash)"
    micromamba activate python3.12
else
    echo "Neither conda nor micromamba found in PATH."
    exit 1
fi

# Execute the two commands in parallel
scripts/bedMethyl2cov.py \
    -i output_methylong-5mCG_5hmCG-traditional/ont/ \
    -o output_methylong-5mCG_5hmCG-traditional/cov/ \
    &> 5mCG_5hmCG-traditional-cov.log

sleep 5

scripts/bedMethyl2cov.py \
    -i output_methylong-5mCG_5hmCG-traditional/ont/ \
    -o output_methylong-5mCG_5hmCG-traditional/cov_5X/ --custom_score 5 \
    &> 5mCG_5hmCG-traditional-cov_5X.log
