
library(readr)
library(dplyr)
library(IRanges)
library(GenomicRanges)

read_bedmethyl <- function(
    bedmethyl_file,
    sample,
    columns=c("valid_coverage", "percent_modified"),
    filter_regions = NULL,
    n_max = Inf
  ) {
  col_names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "color", "valid_coverage", "percent_modified", "count_modified", "count_canonical", "count_other_mod", "count_delete", "count_fail", "count_diff", "count_nocall")
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

  if (!is.null(filter_regions)) {
    hits <- GenomicRanges::findOverlaps(gr_methylation, filter_regions)
    gr_methylation <- gr_methylation[S4Vectors::queryHits(hits)]
  }

  results <- list(
    sample = sample,
    gr_methylation = gr_methylation,
    bedmethyl_file = bedmethyl_file,
    nmax = n_max
  )

  return(results)
}