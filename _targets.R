# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)
library(parallelly)
library(here)
library(readr)
library(data.table)
library(rtracklayer)

# Set target options:
tar_option_set(
  # Packages that your targets need for their tasks.
  packages = c(
    "readr",
    "dplyr",
    "IRanges",
    "GenomicRanges",
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
    name = samplesheet,
    command = {
      metadata <- readr::read_csv(here::here("conf/samplesheet-wf-basecalling.csv"))
      df <- data.frame(
        sample = metadata$alias,
        breed = ifelse(
          grepl("^N", metadata$alias), "Nellore",
          ifelse(grepl("^A", metadata$alias), "Angus", NA)
        ),
        path = here(
          "output_methylong-5mC_5hmC-traditional",
          "ont",
          metadata$alias,
          "pileup",
          paste0(metadata$alias, ".bed.gz")
        )
      )
      df
    }
  ),
  tar_target(
    name = gpg_buffer,
    command = rtracklayer::import.bed(
      here::here("bed/ARS-UCD2.0_CpG-Islands_unix_buffer_merged.bed")
    )
  ),
  tar_target(
    name = bedmethyl_file,
    command = list(
      sample = samplesheet$sample,
      path = samplesheet$path
    ),
    pattern = map(samplesheet)
  ),
  tar_target(
    name = bedmethyl_list,
    command = {
      bedmethyl <- read_bedmethyl(
        bedmethyl_file = bedmethyl_file$path,
        sample = bedmethyl_file$sample,
        filter_regions = gpg_buffer
        # n_max = 1000 # debug
      )
    },
    pattern = map(bedmethyl_file),
    iteration = "list"
  ),
  tar_target(
    name = coverage_data,
    command = {
      # Combine the results from all bedmethyl targets into a single data frame
      df_list <- lapply(bedmethyl_list, function(x) {
        dt <- data.table::data.table(
          sample = x$sample,
          valid_coverage = x$gr_methylation$valid_coverage,
          percent_modified = x$gr_methylation$percent_modified
        )

        # Sample 10,000 rows from each data.table
        dt <- dt[sample(.N, 10000)]

        # Add breed information based on sample name
        dt[, breed := ifelse(
          grepl("^N", sample), "Nellore",
          ifelse(grepl("^A", sample), "Angus", NA)
        )]
      })
      # Combine the list of data.tables into a single data.table
      data.table::rbindlist(df_list, use.names = TRUE, fill = TRUE)
    }
  ),
  tar_quarto(
    name = cpg_traditional_report,
    path = "analysis/03-CpG-traditional.qmd",
    quiet = FALSE
  )
)
