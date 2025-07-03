
source("renv/activate.R")

library(duckdb)
library(DBI)
library(here)

# this is the path to the file
file_path <- here::here(
  "output_methylong_combined/ont/A19_jun/pileup",
  "A19_jun.bed.gz"
)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

# define a query to create the table: exclude columns not needed
query_create <- sprintf(
  "CREATE TABLE methylation AS
   SELECT
     \"#chrom\",
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
  dbQuoteString(con, file_path)
)

# Crea tabella caricando dati
dbExecute(con, query_create)

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
  FROM methylation
")

cat("valid_coverage\n")
cat("   Min. :", summary_valid_coverage$min_valid_coverage, "\n")
cat("   1st Qu.:", summary_valid_coverage$first_quartile, "\n")
cat("   Median :", summary_valid_coverage$median_valid_coverage, "\n")
cat("   Mean   :", summary_valid_coverage$mean_valid_coverage, "\n")
cat("   3rd Qu.:", summary_valid_coverage$third_quartile, "\n")
cat("   Max.   :", summary_valid_coverage$max_valid_coverage, "\n")
cat("   SD     :", summary_valid_coverage$sd_valid_coverage, "\n")
cat("   N      :", summary_valid_coverage$n, "\n")

dbDisconnect(con, shutdown=TRUE)
