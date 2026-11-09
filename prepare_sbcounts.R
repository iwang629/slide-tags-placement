# Prepare cell–bead links from a Macosko-Pipelines Slide-tags SBcounts.h5 file.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript prepare_sbcounts.R SBcounts.h5 called_cells.txt output_dir [neighbor_index=36]\n")
  quit(status = 0)
}
if (!length(args) %in% 3:4) stop("See --help")
suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("rhdf5", quietly = TRUE) || !requireNamespace("RANN", quietly = TRUE))
  stop("Install rhdf5 and RANN to read SBcounts.h5")
sb_file <- args[1L]
cell_file <- args[2L]
out_dir <- args[3L]
knn <- if (length(args) == 4L) as.integer(args[4L]) else 36L
stopifnot(file.exists(sb_file), file.exists(cell_file), is.finite(knn), knn >= 2L)
if (file.exists(out_dir)) stop("Output path already exists: ", out_dir)
fetch <- function(path) rhdf5::h5read(sb_file, path)
cb_list <- as.character(fetch("lists/cb_list"))
sb_list <- as.character(fetch("lists/sb_list"))
dt <- data.table(
  cb_fuzzy = cb_list[as.integer(fetch("matrix/cb_index"))],
  umi_id = as.integer(fetch("matrix/umi")),
  sb = sb_list[as.integer(fetch("matrix/sb_index"))],
  reads = as.integer(fetch("matrix/reads"))
)
stopifnot(!anyNA(dt), all(dt$reads > 0L))

cells_original <- readLines(cell_file, warn = FALSE)
cells_original <- sub("-[0-9]+$", "", cells_original)
stopifnot(length(cells_original) > 0L, !anyDuplicated(cells_original),
  uniqueN(nchar(cells_original)) == 1L,
  all(grepl("^[ACGT]+$", cells_original)))

remap_10x <- function(x) {
  basemap <- c(TC = "AG", AG = "TC", GT = "CA", CA = "GT")
  paste0(substr(x, 1L, 7L), unname(basemap[substr(x, 8L, 9L)]),
    substr(x, 10L, 16L))
}
can_remap <- all(nchar(cells_original) == 16L) &&
  all(unique(substr(cells_original, 8L, 9L)) %in% c("TC", "AG", "GT", "CA"))
reads_noremap <- dt[cb_fuzzy %chin% cells_original, sum(reads)]
cells_remapped <- if (can_remap) remap_10x(cells_original) else cells_original
reads_remap <- if (can_remap) dt[cb_fuzzy %chin% cells_remapped, sum(reads)] else -Inf
remap <- can_remap && reads_remap > reads_noremap
cells_match <- if (remap) cells_remapped else cells_original

# Reproduce the upstream HD1 matching rule: exact called barcodes win; a
# non-whitelist observed barcode is retained only when it has one called-cell
# neighbor at Hamming distance one. N is included as an alternate base because
# that is the behavior of the official helper.
observed <- unique(dt$cb_fuzzy)
bases <- c("A", "C", "G", "T", "N")
hd1_parts <- vector("list", length(cells_match))
for (j in seq_along(cells_match)) {
  cb <- cells_match[j]
  chars <- strsplit(cb, "", fixed = TRUE)[[1L]]
  z <- unlist(lapply(seq_along(chars), function(pos) {
    alts <- bases[bases != chars[pos]]
    paste0(substr(cb, 1L, pos - 1L), alts, substr(cb, pos + 1L, nchar(cb)))
  }), use.names = FALSE)
  z <- z[!z %chin% cells_match & z %chin% observed]
  if (length(z)) hd1_parts[[j]] <- data.table(cb_fuzzy = z, cb_match = cb)
}
fuzzy_map <- rbindlist(hd1_parts, use.names = TRUE, fill = TRUE)
if (nrow(fuzzy_map)) {
  fuzzy_map[, mapping_n := uniqueN(cb_match), by = cb_fuzzy]
  fuzzy_map <- unique(fuzzy_map[mapping_n == 1L, .(cb_fuzzy, cb_match)])
}
exact_map <- data.table(cb_fuzzy = cells_match[cells_match %chin% observed],
  cb_match = cells_match[cells_match %chin% observed])
cb_map <- rbindlist(list(fuzzy_map, exact_map), use.names = TRUE, fill = TRUE)
stopifnot(!anyDuplicated(cb_map$cb_fuzzy))
dt <- merge(dt, cb_map, by = "cb_fuzzy", all.x = TRUE, sort = FALSE)
dt <- dt[!is.na(cb_match)]
dt[, cb := if (remap) remap_10x(cb_match) else cb_match]
stopifnot(dt$cb %chin% cells_original)

# Reproduce ReadPuck filtering and multiplexed-puck x stacking.
puck <- data.table(
  sb = as.character(fetch("puck/sb")),
  x = as.numeric(fetch("puck/x")),
  y = as.numeric(fetch("puck/y")),
  puck_index = as.integer(fetch("puck/puck_index"))
)
stopifnot(!anyNA(puck))
chars <- strsplit(puck$sb, "", fixed = TRUE)
puck[, most_character := vapply(chars, function(v) max(tabulate(match(v, c("A","C","G","T","N")), 5L)), integer(1))]
puck[, longest_run := vapply(chars, function(v) {
  if (!length(v)) return(0L)
  max(rle(v)$lengths)
}, integer(1))]
sb_len <- max(nchar(puck$sb))
mc_tol <- round(sb_len * 4 / 5)
lr_tol <- round(sb_len * 5 / 7)
dups <- unique(puck[duplicated(sb) | duplicated(sb, fromLast = TRUE), sb])
puck <- puck[!sb %chin% dups & !grepl("N", sb) &
  most_character <= mc_tol & longest_run <= lr_tol]
puck[, c("most_character", "longest_run") := NULL]
puck_parts <- split(puck, by = "puck_index", sorted = TRUE, keep.by = FALSE)
if (length(puck_parts) > 1L) {
  widths <- vapply(puck_parts, function(z) max(z$x) - min(z$x), numeric(1))
  starts <- c(0, head(cumsum(widths), -1L))
  puck_parts <- Map(function(z, start) {
    z[, x := x - min(x) + start]
    z
  }, puck_parts, starts)
}
puck <- rbindlist(puck_parts)
stopifnot(!anyDuplicated(puck$sb), nrow(puck) >= knn)

# Match the official post-CB processing: retain valid puck beads, collapse
# reads per cell/UMI/bead, keep a UMI only when one bead has a unique maximum
# read count, and count retained UMI identities per cell/bead.
dt <- dt[sb %chin% puck$sb]
dt <- dt[, .(reads = sum(reads)), by = .(cb, umi_id, sb)]
dt[, keep := reads == max(reads) & sum(reads == max(reads)) == 1L,
  by = .(cb, umi_id)]
dt <- dt[keep == TRUE, .(umi = .N), by = .(cb, sb)]
links <- merge(dt, puck, by = "sb", all = FALSE, sort = FALSE)
setcolorder(links, c("cb", "umi", "sb", "x", "y"))
setorder(links, cb, -umi, sb)

set.seed(1)
query_rows <- sample(seq_len(nrow(puck)), min(10000L, nrow(puck)))
eps <- median(RANN::nn2(data = puck[, .(x, y)], query = puck[query_rows, .(x, y)],
  k = knn)$nn.dists[, knn])
stopifnot(is.finite(eps), eps > 0)

dir.create(out_dir, recursive = TRUE)
fwrite(links[, .(cb, sb, umi, x, y)], file.path(out_dir, "links.csv.gz"))
fwrite(puck[, .(sb, x, y)], file.path(out_dir, "puck_beads.csv.gz"))
writeLines(cells_original, file.path(out_dir, "cells.txt"))
fwrite(data.table(eps_um = eps, neighbor_index = knn, sampled_beads = length(query_rows),
  called_cells = length(cells_original), remapped = remap, puck_beads = nrow(puck),
  links = nrow(links)), file.path(out_dir, "metadata.csv"))
