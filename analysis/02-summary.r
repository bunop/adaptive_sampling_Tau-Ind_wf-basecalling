
source("renv/activate.R")

library(duckdb)
library(DBI)
library(here)

sample_dir <- "output_methylong-5mC_5hmC/ont/A19_jun/pileup"

# this is the path to the file
methylation_path <- here::here(
  sample_dir,
  "A19_jun.bed.gz"
)
regions_path <- here::here(
  "bed",
  "ARS-UCD2.0_CpG-Islands_unix_buffer_merged.bed"
)

# attempt to save data to a file
db_dir <- here::here("analysis", "A19_jun.duckdb")
con <- dbConnect(duckdb::duckdb(), dbdir = db_dir)

# create the table if it does not exist
if (!dbExistsTable(con, "methylation")) {
  create_methylation <- sprintf(
    "CREATE TABLE methylation AS
     SELECT
       \"#chrom\" AS chrom,
       chromStart,
       chromEnd,
       name,
       strand,
       valid_coverage,
       percent_modified,
       count_modified,
       count_canonical,
       count_other_mod,
       count_delete,
       count_fail,
       count_diff,
       count_nocall
     FROM read_csv_auto(%s, delim='\t', header=TRUE)",
    dbQuoteString(con, methylation_path)
  )
  dbExecute(con, create_methylation)
  cat("Table 'methylation' created successfully.\n")
}

# create a new table to load bed file
if (!dbExistsTable(con, "regions")) {
  create_regions <- sprintf(
    "CREATE TABLE regions AS
     SELECT
       column0 AS chrom,
       column1 AS chromStart,
       column2 AS chromEnd,
     FROM read_csv_auto(%s, delim='\t', header=FALSE)",
    dbQuoteString(con, regions_path)
  )
  dbExecute(con, create_regions)
  cat("Table 'regions' created successfully.\n")
}

# crate a table with methylation data in regions
if (!dbExistsTable(con, "methylation_regions")) {
  create_methylation_regions <- "
    CREATE TABLE methylation_regions AS
      SELECT
        m.*
      FROM methylation AS m
      INNER JOIN regions AS r
      ON m.chrom = r.chrom
        AND m.chromStart >= r.chromStart
        AND m.chromEnd < r.chromEnd
  "
  dbExecute(con, create_methylation_regions)
  cat("Table 'methylation_regions' created successfully.\n")
}

# define a summary of a column
summary_valid_coverage <- dbGetQuery(con, "
  SELECT
    COUNT(*) AS n,
    MIN(valid_coverage) AS min_valid_coverage,
    MAX(valid_coverage) AS max_valid_coverage,
    AVG(valid_coverage) AS mean_valid_coverage,
    MEDIAN(valid_coverage) AS median_valid_coverage,
    STDDEV(valid_coverage) AS sd_valid_coverage,
    QUANTILE_CONT(valid_coverage, 0.25) AS first_quartile,
    QUANTILE_CONT(valid_coverage, 0.75) AS third_quartile
  FROM methylation_regions
  USING SAMPLE 1 PERCENT (bernoulli)
")

cat("valid_coverage\n")
cat("   Min.   :", summary_valid_coverage$min_valid_coverage, "\n")
cat("   1st Qu.:", summary_valid_coverage$first_quartile, "\n")
cat("   Median :", summary_valid_coverage$median_valid_coverage, "\n")
cat("   Mean   :", summary_valid_coverage$mean_valid_coverage, "\n")
cat("   3rd Qu.:", summary_valid_coverage$third_quartile, "\n")
cat("   Max.   :", summary_valid_coverage$max_valid_coverage, "\n")
cat("   SD     :", summary_valid_coverage$sd_valid_coverage, "\n")
cat("   N      :", summary_valid_coverage$n, "\n")

dbDisconnect(con, shutdown = TRUE)
