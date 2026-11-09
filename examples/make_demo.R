#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
out <- commandArgs(trailingOnly = TRUE)
if (length(out) != 1L || file.exists(out)) stop("Usage: Rscript examples/make_demo.R new_demo_dir")
dir.create(out, recursive = TRUE)
manifest <- rbindlist(lapply(seq_len(2L), function(k) {
  run <- paste0("lane", k)
  cells <- sprintf("cell%03d", seq_len(28L))
  rna <- data.table(cb = cells, rna_umi = 100 + seq_along(cells) * 80 + k)
  # Reusing a barcode in a different lane does not identify the same cell.
  links <- rbindlist(lapply(seq_len(26L), function(i) {
    x <- 140 * i
    a <- data.table(cb = cells[i], sb = paste0("p", i, "-", 1:4),
      umi = rep(if (i == 26L) 7L else 11L + i, 4),
      x = x + c(0, 2, 0, 2), y = c(0, 0, 2, 2))
    b <- data.table(cb = cells[i], sb = paste0("q", i, "-", 1:3), umi = c(6L, 7L, 7L),
      x = x + 35 + c(0, 2, 0), y = c(0, 0, 2))
    rbind(a, b, data.table(cb = cells[i], sb = paste0("weak", i),
      umi = 1L, x = x + 65, y = 0))
  }))
  # Include a separate, highly dense patch with links too weak to place a cell.
  dense <- CJ(x = 1:6, y = 1:6)
  dense[, `:=`(cb = cells[27], sb = paste0("dense", seq_len(.N)),
    umi = 1L, x = x / 2, y = y / 2)]
  fwrite(rbind(links, dense, use.names = TRUE), file.path(out, paste0(run, "_links.csv.gz")))
  fwrite(rna, file.path(out, paste0(run, "_rna.csv")))
  writeLines(cells, file.path(out, paste0(run, "_cells.txt")))
  data.table(sample = "demo-puck", run = run, links = paste0(run, "_links.csv.gz"),
    cells = paste0(run, "_cells.txt"), rna = paste0(run, "_rna.csv"), eps_um = 8)
}))
fwrite(manifest, file.path(out, "manifest.csv"))
