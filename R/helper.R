
library(readr)
library(dplyr)
library(stringr)
library(IRanges)
library(GenomicRanges)
library(data.table)

# define a test region for simplicity
TEST_REGION <- GenomicRanges::GRanges(
  seqname = "NC_037328.1",
  ranges = IRanges::IRanges(start = 1, end = 2 * 10^7)
)

read_samplesheet <- function(samplesheet_file, results_dir) {
  metadata <- readr::read_csv(here::here(samplesheet_file))

  df <- data.frame(
    sample = metadata$alias,
    breed = ifelse(
      grepl("^N", metadata$alias), "Nellore",
      ifelse(grepl("^A", metadata$alias), "Angus", NA)
    ),
    path = here(
      results_dir,
      "ont",
      metadata$alias,
      "pileup",
      paste0(metadata$alias, ".bed.gz")
    )
  )

  return(df)
}

read_bedmethyl <- function(
    bedmethyl_file,
    sample,
    columns = c("name", "valid_coverage", "percent_modified"),
    filter_regions = NULL,
    n_max = Inf,
    debug = FALSE
  ) {

  # define column names and types for reading the BED file
  col_names <- c(
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "color",
    "valid_coverage",
    "percent_modified",
    "count_modified",
    "count_canonical",
    "count_other_mod",
    "count_delete",
    "count_fail",
    "count_diff",
    "count_nocall"
  )

  col_types <- readr::cols(
    chrom = readr::col_character(),
    chromStart = readr::col_integer(),
    chromEnd = readr::col_integer(),
    name = readr::col_character(),
    score = readr::col_integer(),
    strand = readr::col_character(),
    thickStart = readr::col_integer(),
    thickEnd = readr::col_integer(),
    color = readr::col_character(), # RGB triplet as string
    valid_coverage = readr::col_integer(),
    percent_modified = readr::col_double(),
    count_modified = readr::col_integer(),
    count_canonical = readr::col_integer(),
    count_other_mod = readr::col_integer(),
    count_delete = readr::col_integer(),
    count_fail = readr::col_integer(),
    count_diff = readr::col_integer(),
    count_nocall = readr::col_integer()
  )

  # Read the BED file with methylation data
  bedMethyl <- readr::read_tsv(bedmethyl_file, col_names = col_names, col_types = col_types, comment = "#", n_max = n_max) %>%
    dplyr::select(chrom, chromStart, chromEnd, strand, dplyr::all_of(columns)) %>%
    dplyr::mutate(strand = ifelse(strand == ".", "*", strand))

  gr_methylation <- GenomicRanges::GRanges(
    seqnames = bedMethyl$chrom,
    ranges = IRanges::IRanges(start = bedMethyl$chromStart + 1, end = bedMethyl$chromEnd),
    strand = bedMethyl$strand
  )

  # Add metadata columns
  mcols(gr_methylation) <- bedMethyl %>%
    dplyr::select(dplyr::all_of(columns))

  # debug: take a test region
  if (debug) {
    gr_methylation <- IRanges::subsetByOverlaps(
      gr_methylation,
      TEST_REGION
    )
  }

  # if I provide a Granges object with regions to filter,
  # subset the methylation data
  if (!is.null(filter_regions)) {
    gr_methylation <- IRanges::subsetByOverlaps(
      gr_methylation,
      filter_regions
    )
  }

  results <- list(
    sample = sample,
    gr_methylation = gr_methylation,
    bedmethyl_file = bedmethyl_file,
    nmax = n_max
  )

  return(results)
}

# Function to summarize coverage data from a list of bedmethyl objects
# It extracts valid coverage and percent modified from each bedmethyl object
# and returns a summary data.table with these statistics.
# The function assumes that each bedmethyl object has a gr_methylation slot
# containing valid_coverage and percent_modified fields.
summarize_coverage_data <- function(bedmethyl_list) {
  df_list <- lapply(bedmethyl_list, function(x) {
    dt <- data.table::data.table(
      valid_coverage = x$gr_methylation$valid_coverage,
      percent_modified = x$gr_methylation$percent_modified
    )
  })

  # Combine the list of data.tables into a single data.table
  dt <- data.table::rbindlist(df_list, use.names = TRUE, fill = TRUE)

  # Calculate summary statistics
  return(summary(dt))
}

get_coverage_data <- function(bedmethyl_list, model, n_subsample = NULL) {
  df_list <- lapply(bedmethyl_list, function(x) {
    dt <- data.table::data.table(
      sample = x$sample,
      name = x$gr_methylation$name,
      valid_coverage = x$gr_methylation$valid_coverage,
      percent_modified = x$gr_methylation$percent_modified
    )

    # Subsampling balanced by name
    # TODO: take the same amount of rows per name?
    if (!is.null(n_subsample) && nrow(dt) > n_subsample) {
      # Count the number of unique names
      n_name <- dt[, .N, by = name]

      # Calculate the number of rows to sample per name
      # ie. two names -> n_subsamples / 2 elements for each name
      n_per_name <- floor(n_subsample / nrow(n_name))

      # Sample rows per name
      dt <- dt[, .SD[sample(.N, min(.N, n_per_name))], by = name]
    }

    # Add breed information based on sample name
    dt[, breed := ifelse(
      grepl("^N", sample), "Nellore",
      ifelse(grepl("^A", sample), "Angus", NA)
    )]
  })

  # Combine the list of data.tables into a single data.table
  dt <- data.table::rbindlist(df_list, use.names = TRUE, fill = TRUE)

  # Add model label
  dt$model <- model

  return(dt)
}

load_bsseq <- function(samplesheet, cpg_buffer, debug = FALSE) {
  # define metadata for BSseq
  metadata <- samplesheet %>%
    dplyr::select(sample, breed) %>%
    tibble::column_to_rownames("sample")

  # read the bedmethyl files (no motif as names)
  BS.seq <- bsseq::read.modkit(
    samplesheet$path,
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
  if (debug) {
    BS.seq <- IRanges::subsetByOverlaps(
      BS.seq,
      GenomicRanges::GRanges(
        seqname = "NC_037328.1",
        ranges = IRanges::IRanges(start = 1, end = 2 * 10^7)
      )
    )
  }

  return(BS.seq)
}

bsseq_genetrack <- function(gff_exons) {
  # Convert GFF exons to a data frame for BSseq
  df <- data.frame(
    chr = as.character(seqnames(gff_exons)),
    start = start(gff_exons),
    end = end(gff_exons),
    strand = as.character(strand(gff_exons)),
    gene_ID = gff_exons$ID,
    gene_name = gff_exons$gene,
    transcript_id = gff_exons$transcript_id
  )

  # counting isoforms and exon numbers
  df <- df %>%
    dplyr::mutate(
      # remove versioning from transcript_id
      transcript_base = stringr::str_remove(transcript_id, "-\\d+$")
    ) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(
      exon_number = dplyr::row_number() - 1,                 # start from 0
      isoforms = as.integer(factor(transcript_base))  # progressive isoform
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      isoforms = ifelse(is.na(isoforms), 1L, isoforms)
    )

  # remove unnecessary columns
  df <- df %>% select(-transcript_id, -transcript_base)

  return(df)
}