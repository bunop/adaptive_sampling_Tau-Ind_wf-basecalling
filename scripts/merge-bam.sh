#!/bin/bash
#SBATCH --nodes=1                                               # 1 node
#SBATCH --ntasks-per-node=1                                     # 1 tasks per node
#SBATCH --cpus-per-task=4                                       # x CPUs per task
#SBATCH --time=24:00:00                                         # time limits: see queue and QoS
#SBATCH --mem=16G                                               # 16GB to manage process
#SBATCH --error=merge-bam.err                                   # standard error file
#SBATCH --output=merge-bam.out                                  # standard output file
#SBATCH --job-name=merge-bam                                    # job name
#SBATCH --account=IscrC_BioGPUPX                                # account name
#SBATCH --partition=g100_usr_prod                               # partition name (see https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.3%3A+GALILEO100+UserGuide)
#SBATCH --qos=normal                                            # quality of service (see https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.3%3A+GALILEO100+UserGuide)

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
