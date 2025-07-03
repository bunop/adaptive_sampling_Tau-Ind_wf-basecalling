#!/bin/bash
#SBATCH --nodes=1                                               # 1 node
#SBATCH --ntasks-per-node=1                                     # 1 tasks per node
#SBATCH --cpus-per-task=1                                       # 2 CPUs per task
#SBATCH --time=4-00:00:00                                       # time limits: see queue and QoS
#SBATCH --mem=16G                                               # 16GB to manage process
#SBATCH --error=methylong.err                                   # standard error file
#SBATCH --output=methylong.out                                  # standard output file
#SBATCH --job-name=methylong                                    # job name
#SBATCH --account=IscrC_BioGPUPX                                # account name
#SBATCH --partition=g100_usr_prod                               # partition name (see https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.3%3A+GALILEO100+UserGuide)
#SBATCH --qos=g100_qos_lprod                                    # quality of service (see https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.3%3A+GALILEO100+UserGuide)

# set the path of institution-specific configuration files
# (required since we're working offline)
export CUSTOM_CONFIG_BASE=${WORK}/nf-configs

# mind to the pipeline version (required)
nextflow run nf-core/methylong -r 1.0.0 \
    --custom_config_base ${CUSTOM_CONFIG_BASE} \
    -config ${CUSTOM_CONFIG_BASE}/nfcore_custom.config \
    -config conf/custom-methylong.config \
    -profile ibba,galileo -resume -params-file conf/params-methylong.json
