
# TODO

- [ ] load coverage data from `bsseq` objects
- [x] display summaries on all *bedMethyl files*
- [ ] plot coverage and % modified *after* filtering for 5X coverage
- [ ] display `bsseq` object in report (the number of loci selected)
- [ ] describe `dmrs` objects
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

- [ ] deal with `HDF5` objects
