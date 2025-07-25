
# TODO

- [ ] load coverage data from `bsseq` objects
- [ ] Plot a hexbin plot of percent modified methylation vs coverage:

  ```{r plot-percent-modified-coverage}
  ggplot(coverage_data, aes(x = score, y = percent_modified)) +
    geom_hex(bins = 50) +
    scale_x_log10() +
    scale_fill_viridis_c(trans = "log10") +  # <--- Scala logaritmica sui conteggi
    labs(
      title = "Hexbin: Percent Modified vs Coverage (log10 scale)",
      x = "Coverage (log10 scale)",
      y = "Methylation (%)",
      fill = "Count (log10)"
    ) +
    theme_minimal()
  ```
