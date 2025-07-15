#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --error=merge-bam-5mC_5hmC.err
#SBATCH --output=merge-bam-5mC_5hmC.out
#SBATCH --job-name=merge-bam-5mC_5hmC
#SBATCH --account=IscrC_BioGPUPX    # account name
#SBATCH --partition=g100_usr_prod   # partition name (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)
#SBATCH --qos=normal                # quality of service (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)

set -euo pipefail

SAMTOOLS="singularity run -B ${HOME} -B ${WORK} -B ${PWD} -B ${SCRATCH} ${NXF_SINGULARITY_CACHEDIR}/samtools-1.21--h50ea8bc_0.img samtools"
CPUS=${SLURM_CPUS_PER_TASK:-1}

for pattern in "A19_jun" "A21_jun" "A25_jun" "N03_jun" "N07_jun" "N13_jun" ; do
    if [ ! -f data/samples-5mC_5hmC/${pattern}.bam ]; then
        $SAMTOOLS merge -@ ${CPUS} -o data/samples-5mC_5hmC/${pattern}.bam demux-5mC_5hmC/*${pattern}.bam
        echo "data/samples-5mC_5hmC/${pattern}.bam created"
    else
        echo "Skipping data/samples-5mC_5hmC/${pattern}.bam (already exists)"
    fi
done

echo "All BAM files have been processed. Program finished."
