
# Adaptive sampling Tau-Ind

Basecalling workflow using `wf-basecalling` and `nf-core/methylong`
pipelines for the `adaptive_sampling_Tau-Ind` dataset. Using
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` model for basecalling.
We will do a basecalling using different models, we will have multiple
configuration files for each context.

## Update models

To update the models, run the following command:

```bash
singularity run ${NXF_SINGULARITY_CACHEDIR}/ontresearch-dorado-shae9327ad17e023b76e4d27cf287b6b9d3a271092b.img \
    dorado download --models-directory ${HOME}/Projects/dorado_models/
```

Then ensure that the `conf/custom-wf-basecalling.config` file has the correct path to the models:

```groovy
env {
    DRD_MODELS_PATH = '${HOME}/Projects/dorado_models'
}

singularity {
    runOptions = '-B ${HOME}/Projects/dorado_models'
}
```

## Calling wf-basecalling pipeline

Called the latest `wf-basecalling` pipeline (`v1.5.2`) customize to support
different models with different context (https://github.com/bunop/wf-basecalling/tree/multiple_calling)
using:

```bash
nextflow run ~/Projects/wf-basecalling/ -profile singularity -resume \
    -c conf/custom-wf-basecalling.config -params-file conf/params-wf-basecalling-5mC_5hmC.json
nextflow run ~/Projects/wf-basecalling/ -profile singularity -resume \
    -c conf/custom-wf-basecalling.config -params-file conf/params-wf-basecalling-5mCG_5hmCG.json
```

> NOTE: is not possible to call `--duplex=true` and `--barcode_kit=SQK-NBD114-24`
> at the same time, *demultiplexing* should be done after basecalling.

## Post processing

### calling PyCOQC

Calculate quality control metrics using `pycoQC`:

```bash
singularity run $NXF_SINGULARITY_CACHEDIR/pip_pycoqc_setuptools_31d5a8754dcc1b68.sif \
    pycoQC -f output-5mC_5hmC/SAMPLE.summary.tsv.gz -o output-5mC_5hmC/SAMPLE.summary.html
singularity run $NXF_SINGULARITY_CACHEDIR/pip_pycoqc_setuptools_31d5a8754dcc1b68.sif \
    pycoQC -f output-5mCG_5hmCG/SAMPLE.summary.tsv.gz -o output-5mCG_5hmCG/SAMPLE.summary.html
```

### Join passed simplex and duplex reads

Called reads are divided into two groups: `simplex` and `duplex` for both passed
and failed reads:

```text
SAMPLE.fail.duplex.cram
SAMPLE.fail.simplex.cram
SAMPLE.pass.duplex.cram
SAMPLE.pass.simplex.cram
```

Let's join the `pass` reads into a single file:

```bash
cd output-5mC_5hmC
singularity run $NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img \
    samtools merge --threads 4 -o SAMPLE.pass.all.cram SAMPLE.pass.duplex.cram SAMPLE.pass.simplex.cram
singularity run $NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img \
    samtools index --threads 4 SAMPLE.pass.all.cram
cd ..
cd output-5mCG_5hmCG
singularity run $NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img \
    samtools merge --threads 4 -o SAMPLE.pass.all.cram SAMPLE.pass.duplex.cram SAMPLE.pass.simplex.cram
singularity run $NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img \
    samtools index --threads 4 SAMPLE.pass.all.cram
cd ..
```

### Demultiplexing

Do demultiplexing using `dorado`:

```bash
singularity run $NXF_SINGULARITY_CACHEDIR/ontresearch-dorado-shae9327ad17e023b76e4d27cf287b6b9d3a271092b.img \
    dorado demux --kit-name SQK-NBD114-24 --threads 4 --verbose --output-dir demux-5mC_5hmC \
    --sample-sheet conf/samplesheet-wf-basecalling.csv output-5mC_5hmC/SAMPLE.pass.all.cram
singularity run $NXF_SINGULARITY_CACHEDIR/ontresearch-dorado-shae9327ad17e023b76e4d27cf287b6b9d3a271092b.img \
    dorado demux --kit-name SQK-NBD114-24 --threads 4 --verbose --output-dir demux-5mCG_5hmCG \
    --sample-sheet conf/samplesheet-wf-basecalling.csv output-5mCG_5hmCG/SAMPLE.pass.all.cram
```

## Call nf-core/methylong

### Preparing samples

The three different runs are in the `demux-5mC_5hmC` directory, so we need to prepare
the samples for `nf-core/methylong`. Let's create a new directory:

```bash
mkdir -p $SCRATCH/adaptive_sampling_Tau-Ind_wf-basecalling/samples-5mC_5hmC/
ln -s $SCRATCH/adaptive_sampling_Tau-Ind_wf-basecalling/samples-5mC_5hmC data/
sbatch scripts/merge-bam-5mC_5hmC.sh
mkdir -p $SCRATCH/adaptive_sampling_Tau-Ind_wf-basecalling/samples-5mCG_5hmCG/
ln -s $SCRATCH/adaptive_sampling_Tau-Ind_wf-basecalling/samples-5mCG_5hmCG data/
sbatch scripts/merge-bam-5mCG_5hmCG.sh
```

### Running nf-core/methylong

Run the `nf-core/methylong` pipeline using the following command:

```bash
sbatch scripts/launch-methylong-5mC_5hmC.sh
sbatch scripts/launch-methylong-5mC_5hmC-cpg.sh
sbatch scripts/launch-methylong-5mC_5hmC-traditional.sh
sbatch scripts/launch-methylong-5mCG_5hmCG.sh
sbatch scripts/launch-methylong-5mCG_5hmCG-cpg.sh
sbatch scripts/launch-methylong-5mCG_5hmCG-traditional.sh
```
