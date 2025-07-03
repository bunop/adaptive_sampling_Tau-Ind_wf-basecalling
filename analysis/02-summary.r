
source("renv/activate.R")

library(duckdb)
library(DBI)
library(here)

# resto dello script ...

library(duckdb)
library(DBI)
library(here)

file_path <- here::here(
  "output_methylong_combined/ont/A19_jun/pileup",
  "A19_jun.bed.gz"
)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

query_create <- sprintf(
  "CREATE TABLE methylation AS
   SELECT * FROM read_csv_auto(%s, delim='\t', header=TRUE)",
  dbQuoteString(con, file_path)
)

# Crea tabella caricando dati
dbExecute(con, query_create)

# Fai summary della colonna score
summary_score <- dbGetQuery(con, "
  SELECT
    COUNT(*) AS n,
    MIN(score) AS min_score,
    MAX(score) AS max_score,
    AVG(score) AS mean_score,
    MEDIAN(score) AS median_score,
    STDDEV(score) AS sd_score
  FROM methylation
")

print(summary_score)

dbDisconnect(con, shutdown=TRUE)
