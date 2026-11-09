suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("R/clustering.R")
source("R/placement.R")
source("R/layout.R")
args <- commandArgs(trailingOnly = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(..., status = 0L) {
  res <- system2(rscript, shQuote(c(...)), stdout = TRUE, stderr = TRUE)
  actual <- attr(res, "status")
  if (is.null(actual)) actual <- 0L
  if (actual != status) stop(paste(res, collapse = "\n"))
  invisible(res)
}
demo <- tempfile("placement-example-")
out <- tempfile("placement-output-")
run("examples/make_demo.R", demo)
manifest <- file.path(demo, "manifest.csv")
run("run_placement.R", manifest, out, "--cores", "2")
stopifnot(file.exists(file.path(out, "DONE")))
a <- rbindlist(lapply(Sys.glob(file.path(out, "runs", "*_audit.csv.gz")), read_table))
s <- read_table(file.path(out, "sample_summary.csv"))
dm <- read_table(file.path(out, "density", "demo-puck.csv.gz"))
stopifnot(nrow(a) == 56L, !anyDuplicated(a, by = c("run", "cb")),
  all(a[cb == "cell028", status] == "NO_QUALIFYING_CLUSTER"),
  all(a[cb == "cell027", filtered_spatial_umi] == 0),
  all(a[cb == "cell026", status] == "AMBIGUOUS_COMPETITORS"),
  all(a[cb == "cell020", accepted]),
  all(a[cb == "cell020", x] == 2801), all(a[cb == "cell020", y] == 1),
  all(a[accepted == FALSE, is.na(x) & is.na(y)]),
  all(dm[grepl("^dense", sb), density_removed]),
  s$eligible == 56L, s$placed == sum(a$accepted),
  s$adaptive_threshold == 30, s$depth_reference_source == "positive_UMI_fallback")

# Identical saved fold assignments, acceptance and positions on a serial rerun.
folds <- file.path(demo, "folds.csv")
fwrite(a[, .(run, cb, fold)], folds)
again <- tempfile("placement-replay-")
run("run_placement.R", manifest, again, "--folds", folds)
b <- rbindlist(lapply(Sys.glob(file.path(again, "runs", "*_audit.csv.gz")), read_table))
stopifnot(identical(names(a), names(b)), all(vapply(names(a),
  function(nm) identical(a[[nm]], b[[nm]]), logical(1))))

# One lane also works; a puck with no links reports rejection, not a crash.
m <- read_table(manifest)
single <- file.path(demo, "one.csv")
fwrite(m[1L], single)
single_out <- tempfile("placement-single-")
run("run_placement.R", single, single_out)
stopifnot(read_table(file.path(single_out, "sample_summary.csv"))$eligible == 28L)
fwrite(data.table(cb = character(), sb = character(), umi = integer(), x = numeric(), y = numeric()),
  file.path(demo, "empty.csv"))
fwrite(m[1L][, links := "empty.csv"], single)
empty_out <- tempfile("placement-empty-")
run("run_placement.R", single, empty_out)
stopifnot(read_table(file.path(empty_out, "sample_summary.csv"))$placed == 0)

# If RNA depth has no variation, do not silently treat the hard floor as an
# inferred cutoff. The fallback must be explicitly chosen and stays labeled.
sparse <- data.table(sample = "puck", fold = 1:5, rna_umi = 500, qualifying_n = 1L, primary_umi = 50)
stopifnot(all(is.na(rna_thresholds(sparse)$rna_cutoff)),
  all(rna_thresholds(sparse, on_sparse = "floor")$rna_reference_status == "FLOOR_FALLBACK"),
  all(rna_thresholds(sparse, on_sparse = "floor")$rna_cutoff == 10))

# Rigid reflection and translation preserve pairwise distances, and do not
# give an ungrouped cell organized coordinates.
xy <- CJ(x = 1:5, y = 1:5)
p <- xy[, .(run = "test", sample = "puck", cb = paste0("id", seq_len(.N)), x, y)]
g <- p[1:24, .(run, cb, group = "piece1")]
anchors <- p[1:24, .(run, cb, group = "piece1", target_x = -x + 100, target_y = y - 50)]
z <- organize_placements(p, g, anchors)
stopifnot(z$transforms$reflection, z$transforms$rmse_um < 1e-9,
  is.na(z$coordinates$layout_x[25L]), max(abs(z$coordinates$layout_x[1:24] - anchors$target_x)) < 1e-9)
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  counts <- matrix(1L, nrow = 2, ncol = 25, dimnames = list(c("A", "B"), p$cb))
  obj <- SeuratObject::CreateSeuratObject(counts)
  obj2 <- attach_seurat(obj, z$coordinates[1:24], "test")
  stopifnot(identical(colnames(obj2), colnames(obj)), sum(obj2$spatial_placed) == 24L,
    !obj2$spatial_placed[25L], nrow(SeuratObject::Embeddings(obj2[["spatial_layout"]])) == 24L)
}

if (requireNamespace("rhdf5", quietly = TRUE) && requireNamespace("RANN", quietly = TRUE)) {
  h5 <- tempfile(fileext = ".h5")
  rhdf5::h5createFile(h5)
  for (group in c("matrix", "lists", "puck")) rhdf5::h5createGroup(h5, group)
  bs <- c("ACGTACGT", "TGCATGCA", "CGTACGTA", "GTACGTAC")
  cb <- c("ACGT", "TCGT")
  h <- list("lists/cb_list" = cb, "lists/sb_list" = bs,
    "matrix/cb_index" = c(1L, 1L, 1L, 2L, 2L), "matrix/sb_index" = c(1L, 2L, 1L, 3L, 4L),
    "matrix/umi" = c(1L, 1L, 2L, 1L, 1L), "matrix/reads" = c(3L, 1L, 2L, 1L, 1L),
    "puck/sb" = bs, "puck/x" = c(0, 10, 10, 0), "puck/y" = c(0, 0, 10, 10),
    "puck/puck_index" = rep(1L, 4L))
  for (nm in names(h)) rhdf5::h5write(h[[nm]], h5, nm)
  cells <- tempfile()
  writeLines(cb, cells)
  prepared <- tempfile()
  run("prepare_sbcounts.R", h5, cells, prepared, "2")
  l <- read_table(file.path(prepared, "links.csv.gz"))
  stopifnot(nrow(l) == 1L, l$cb == cb[1L], l$sb == bs[1L], l$umi == 2L,
    read_table(file.path(prepared, "metadata.csv"))$eps_um == 10)
}
cat("Synthetic placement, replay, single/empty lane, sparse RNA, and layout checks passed.\n")
