# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)
library(parallelly)

# Set target options:
tar_option_set(
  # Packages that your targets need for their tasks.
  packages = c(
    "here",
    "readr",
    "data.table",
    "dplyr",
    "tibble",
    "bsseq",
    "rtracklayer",
    "IRanges",
    "GenomicRanges",
    "BiocParallel",
    "BiocGenerics",
    "ggplot2",
    "quarto"
  ),
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  controller = crew::crew_controller_local(
    workers = parallelly::availableCores() - 1,
    seconds_idle = 60,
    tasks_max = 50
  ),
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
  seed = 42
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(
    name = samplesheet_5mC_5hmC,
    command = read_samplesheet(
      here::here("conf/samplesheet-wf-basecalling.csv"),
      "output_methylong-5mC_5hmC-traditional"
    )
  ),
  tar_target(
    name = samplesheet_5mCG_5hmCG,
    command = read_samplesheet(
      here::here("conf/samplesheet-wf-basecalling.csv"),
      "output_methylong-5mCG_5hmCG-traditional"
    )
  ),
  # load annotations
  tar_target(
    name = cpg_buffer,
    command = rtracklayer::import.bed(
      here::here("bed/ARS-UCD2.0_CpG-Islands_unix_buffer_merged.bed")
    )
  ),
  tar_target(
    name = geneTrack,
    command = {
      gff_path <- here::here("genome/genomic.gff")
      gff_exons <- rtracklayer::import.gff(gff_path, feature.type = "exon")

      # Convert GFF exons to a data frame for BSseq (helper function)
      bsseq_genetrack(gff_exons)
    }
  ),
  # read bedmethyl files for both models, 5mC_5hmC and 5mCG_5hmCG
  # 5mC_5hmC model
  tar_target(
    name = bedmethyl_file_5mC_5hmC,
    command = list(
      sample = samplesheet_5mC_5hmC$sample,
      path = samplesheet_5mC_5hmC$path
    ),
    pattern = map(samplesheet_5mC_5hmC)
  ),
  tar_target(
    name = bedmethyl_list_5mC_5hmC,
    command = {
      bedmethyl <- read_bedmethyl(
        bedmethyl_file = bedmethyl_file_5mC_5hmC$path,
        sample = bedmethyl_file_5mC_5hmC$sample,
        # n_max = 1000, # debug
        filter_regions = cpg_buffer
      )
    },
    pattern = map(bedmethyl_file_5mC_5hmC),
    iteration = "list"
  ),
  tar_target(
    name = coverage_data_5mC_5hmC,
    command = get_coverage_data(
      bedmethyl_list_5mC_5hmC,
      "5mC_5hmC",
      n_subsample = 100000 # Subsample to 100,000 rows for visualization
    )
  ),
  # 5mCG_5hmCG model
  tar_target(
    name = bedmethyl_file_5mCG_5hmCG,
    command = list(
      sample = samplesheet_5mCG_5hmCG$sample,
      path = samplesheet_5mCG_5hmCG$path
    ),
    pattern = map(samplesheet_5mCG_5hmCG)
  ),
  tar_target(
    name = bedmethyl_list_5mCG_5hmCG,
    command = {
      bedmethyl <- read_bedmethyl(
        bedmethyl_file = bedmethyl_file_5mCG_5hmCG$path,
        sample = bedmethyl_file_5mCG_5hmCG$sample,
        # n_max = 1000, # debug
        filter_regions = cpg_buffer
      )
    },
    pattern = map(bedmethyl_file_5mCG_5hmCG),
    iteration = "list"
  ),
  tar_target(
    name = coverage_data_5mCG_5hmCG,
    command = get_coverage_data(
      bedmethyl_list_5mCG_5hmCG,
      "5mCG_5hmCG",
      n_subsample = 100000 # Subsample to 100,000 rows for visualization
    )
  ),
  # summarize coverage data
  tar_target(
    name = summary_coverage_data_5mC_5hmC,
    command = summary(coverage_data_5mC_5hmC)
  ),
  tar_target(
    name = summary_coverage_data_5mCG_5hmCG,
    command = summary(coverage_data_5mCG_5hmCG)
  ),
  # combine coverage data
  tar_target(
    name = combined_coverage_data,
    command = rbind(coverage_data_5mC_5hmC, coverage_data_5mCG_5hmCG)
  ),
  # Combine coverage data from both models, then make plots
  tar_target(
    name = plot_valid_coverage,
    command = {
      ggplot2::ggplot(combined_coverage_data, aes(y = valid_coverage, x = model, fill = model)) +
        # disable outliers for better visibility
        ggplot2::geom_boxplot(outlier.shape = NA) +
        facet_wrap(~ sample) +
        ggplot2::scale_y_log10(limits = c(NA, 100)) +
        ggplot2::labs(
          title = "Distribution of Valid Coverage in CpG Buffer Regions",
          y = "Valid Coverage (log10 scale)",
        ) +
        ggplot2::theme_minimal()
    }
  ),
  tar_target(
    name = plot_percent_modified,
    command = {
      ggplot2::ggplot(combined_coverage_data, aes(y = percent_modified, x = model, fill = model)) +
        ggplot2::geom_boxplot() +
        facet_wrap(~ sample) +
        ggplot2::labs(
          title = "Distribution of Percent Modified Methylation",
          y = "Modified Methylation (%)"
        ) +
        ggplot2::theme_minimal()
    }
  ),
  # load, sort and filter the BSseq object
  tar_target(
    name = BS.seq,
    command = {
      # define metadata for BSseq
      metadata <- samplesheet_5mC_5hmC %>%
        dplyr::select(sample, breed) %>%
        tibble::column_to_rownames("sample")

      # read the bedmethyl files (no motif as names)
      BS.seq <- bsseq::read.modkit(
        samplesheet_5mC_5hmC$path,
        colData = metadata,
        rmZeroCov = TRUE,
        strandCollapse = TRUE
      )

      # ensure the BSseq object is sorted and filtered by cpg_buffer
      BS.seq.sorted <- bsseq::orderBSseq(BS.seq)
      BS.seq <- IRanges::subsetByOverlaps(
        BS.seq.sorted,
        cpg_buffer
      )

      # debug: take first chromosome
      # BS.seq <- IRanges::subsetByOverlaps(
      #   BS.seq,
      #   GenomicRanges::GRanges(seqname = "NC_037328.1", ranges = IRanges::IRanges(start = 1, end = 2*10^7))
      # )

      return(BS.seq)
    }
  ),
  # smooth the BSseq object
  # TODO: parallelize this step
  # TODO: wait for #149 (https://github.com/hansenlab/bsseq/issues/149) to be resolved
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq_analysis.html#21_Manually_splitting_the_smoothing_computation
  tar_target(
    name = BS.seq.fit,
    command = bsseq::BSmooth(
      BSseq = BS.seq,
      BPPARAM = MulticoreParam(workers = parallelly::availableCores() - 1)
    )
  ),
  # filter low coverage regions
  tar_target(
    name = BS.seq.ex.fit,
    command = {
      BS.cov <- bsseq::getCoverage(BS.seq.fit)

      # keep loci with 5X coverage in at least 2 samples per breed
      keepLoci.ex <- which(
        rowSums(BS.cov[, BS.seq.fit$breed == "Angus"] >= 5) >= 2 &
        rowSums(BS.cov[, BS.seq.fit$breed == "Nellore"] >= 5) >= 2
      )
      BS.seq.fit[keepLoci.ex, ]
    }
  ),
  tar_target(
    name = bsseq_tstat,
    command = {
      # determine the two groups for t-statistic calculation
      metadata <- pData(BS.seq.ex.fit)
      group1 <- rownames(metadata)[metadata$breed == "Angus"]
      group2 <- rownames(metadata)[metadata$breed == "Nellore"]

      bsseq_tstat <- bsseq::BSmooth.tstat(
        BS.seq.ex.fit,
        group1 = group1,
        group2 = group2,
        estimate.var = "same",
        local.correct = TRUE,
        verbose = FALSE
      )

      # need to filter out NA values
      idxs <- which(!is.na(bsseq_tstat@stats[, 6]))
      bsseq_tstat[idxs, ]
    }
  ),
  tar_target(
    name = dmrs,
    command = {
      dmrs0 <- bsseq::dmrFinder(bsseq_tstat)
      BiocGenerics::subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
    }
  ),
  tar_quarto(
    name = cpg_traditional_report,
    path = "analysis/03-CpG-traditional.qmd",
    quiet = FALSE
  )
)
