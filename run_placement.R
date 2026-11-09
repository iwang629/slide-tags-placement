#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || !length(args)) {
  cat(paste0(
    "Usage: Rscript run_placement.R manifest.csv new_output_dir [--option value ...]\n",
    "Required manifest columns: sample,run,links,cells,rna,eps_um (or eps_metadata).\n",
    "Options: --cores 1 --seed 296 --density-quantile 0.90 --link-umi-min 2\n",
    "         --total-floor 30 --total-cap 100 --total-fraction 0.05\n",
    "         --primary-quantile 0.05 --primary-floor 10 --primary-beads 3 --dominance 3\n",
    "         --sparse-rna-reference fail|floor --folds path/to/folds.csv\n"))
  quit(status = 0L)
}
if (length(args) < 2L || (length(args) - 2L) %% 2L != 0L) stop("See --help")
defaults <- c("--cores" = "1", "--seed" = "296", "--density-quantile" = ".90",
  "--link-umi-min" = "2", "--total-floor" = "30", "--total-cap" = "100",
  "--total-fraction" = ".05", "--primary-quantile" = ".05", "--primary-floor" = "10",
  "--primary-beads" = "3", "--dominance" = "3", "--sparse-rna-reference" = "fail", "--folds" = "")
if (length(args) > 2L) {
  keys <- args[seq(3L, length(args), 2L)]
  if (anyDuplicated(keys) || !all(keys %in% names(defaults))) stop("Duplicate or unknown options. See --help")
  defaults[keys] <- args[seq(4L, length(args), 2L)]
}
num <- function(k) {
  x <- suppressWarnings(as.numeric(defaults[[paste0("--", k)]]))
  if (!is.finite(x)) stop("Expected finite number for --", k)
  x
}
cores <- num("cores"); seed <- num("seed")
density_q <- num("density-quantile"); link_min <- num("link-umi-min")
total_floor <- num("total-floor"); total_cap <- num("total-cap"); total_fraction <- num("total-fraction")
primary_q <- num("primary-quantile"); primary_floor <- num("primary-floor")
primary_beads <- num("primary-beads"); dominance <- num("dominance")
sparse <- defaults[["--sparse-rna-reference"]]
stopifnot(cores >= 1, cores == floor(cores), seed >= 0, seed == floor(seed),
  seed <= .Machine$integer.max, density_q > 0, density_q <= 1, link_min >= 1,
  total_floor >= 10, total_cap >= total_floor, total_fraction > 0,
  primary_q > 0, primary_q < 1, primary_floor >= 5,
  primary_beads >= 2, primary_beads == floor(primary_beads), dominance > 1,
  sparse %in% c("fail", "floor"))
suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("dbscan", quietly = TRUE)) stop("Install the R package dbscan")
setDTthreads(1L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
root <- dirname(normalizePath(script))
source(file.path(root, "R", "clustering.R"))
source(file.path(root, "R", "placement.R"))

manifest <- normalizePath(args[1L], mustWork = TRUE)
out <- args[2L]
if (file.exists(out)) stop("Output path already exists: ", out)
m <- read_table(manifest)
require_cols(m, c("sample", "run", "links", "cells", "rna"), "Manifest")
stopifnot(nrow(m) > 0L, valid_ids(m$sample), valid_ids(m$run), !anyDuplicated(m$run),
  all(grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", c(m$sample, m$run))))
path <- function(f) {
  if (is.na(f) || !nzchar(f)) stop("Empty input path")
  normalizePath(if (grepl("^(/|[A-Za-z]:|~)", f)) path.expand(f)
    else file.path(dirname(manifest), f), mustWork = TRUE)
}
for (nm in c("links", "cells", "rna")) m[, (nm) := vapply(get(nm), path, character(1))]
if (!"eps_um" %in% names(m)) m[, eps_um := NA_real_]
if (!"eps_metadata" %in% names(m)) m[, eps_metadata := NA_character_]
for (i in which(is.na(m$eps_um))) {
  e <- read_table(path(m$eps_metadata[i]))
  stopifnot(nrow(e) == 1L, "eps_um" %in% names(e))
  set(m, i, "eps_um", e$eps_um)
}
stopifnot(is.numeric(m$eps_um), all(is.finite(m$eps_um) & m$eps_um > 0))

read_run <- function(i) {
  cells <- readLines(m$cells[i], warn = FALSE)
  if (!length(cells) || !valid_ids(cells) || anyDuplicated(cells) || any(trimws(cells) != cells))
    stop("Expected unique, nonempty barcodes in ", m$cells[i])
  r <- read_table(m$rna[i])
  require_cols(r, c("cb", "rna_umi"), "RNA UMI table")
  if (!valid_ids(r$cb) || anyDuplicated(r$cb) || !is.numeric(r$rna_umi))
    stop("Invalid RNA barcodes or counts: ", m$rna[i])
  rna <- r$rna_umi[match(cells, r$cb)]
  if (any(!is.finite(rna) | rna <= 0 | rna != floor(rna))) stop("Missing or invalid raw RNA UMI for a called cell in ", m$run[i])
  links <- read_table(m$links[i])
  require_cols(links, c("cb", "sb", "umi", "x", "y"), "Link table")
  if (!nrow(links)) links <- data.table(cb = character(), sb = character(),
    umi = numeric(), x = numeric(), y = numeric())
  if (!valid_ids(links$cb) || !valid_ids(links$sb) ||
      !all(vapply(links[, .(umi, x, y)], is.numeric, logical(1))) ||
      any(!is.finite(links$umi) | links$umi <= 0 | links$umi != floor(links$umi)) ||
      any(!is.finite(links$x) | !is.finite(links$y))) stop("Invalid links for ", m$run[i])
  if (!all(links$cb %chin% cells)) stop("Link table has barcodes outside the manifest cell list: ", m$run[i])
  list(links = links[, .(cb, sb, umi, x, y)], cells = cells, rna = rna)
}

fold_override <- NULL
if (nzchar(defaults[["--folds"]])) {
  fold_override <- read_table(defaults[["--folds"]])
  require_cols(fold_override, c("run", "cb", "fold"), "Fold table")
  stopifnot(valid_ids(fold_override$run), valid_ids(fold_override$cb),
    !anyDuplicated(fold_override, by = c("run", "cb")), all(fold_override$run %chin% m$run),
    all(fold_override$fold %in% 1:5))
}
dir.create(out, recursive = TRUE)
dir.create(file.path(out, "density"))
dir.create(file.path(out, "runs"))
summaries <- list()
for (s in unique(m$sample)) {
  message("Processing puck ", s)
  rows <- which(m$sample == s)
  runs <- lapply(rows, read_run)
  bead_parts <- lapply(runs, function(z) unique(z$links[, .(sb, x, y)]))
  pool <- rbindlist(bead_parts)
  dm <- if (nrow(pool)) density_map(pool, median(m$eps_um[rows]), density_q) else
    data.table(sb = character(), x = numeric(), y = numeric(), local_density = integer(),
      density_percentile = numeric(), density_removed = logical())
  fwrite(dm, file.path(out, "density", paste0(s, ".csv.gz")))
  d <- rbindlist(lapply(seq_along(rows), function(j) {
    i <- rows[j]; z <- runs[[j]]
    p <- profile_run(z$links, z$cells, dm, m$eps_um[i], link_min, cores)
    p$audit[, `:=`(sample = s, run = m$run[i], rna_umi = z$rna)]
    p$candidates[, run := rep(m$run[i], .N)]
    fwrite(p$candidates, file.path(out, "runs", paste0(m$run[i], "_candidates.csv.gz")))
    merge(p$audit, summarize_candidates(p$candidates), by = c("run", "cb"), all.x = TRUE, sort = FALSE)
  }))
  d[is.na(qualifying_n), `:=`(qualifying_n = 0L, stable_n = 0L)]
  setorder(d, run, cb)
  if (is.null(fold_override)) {
    set.seed(seed)
    d[, fold := sample(rep(1:5, length.out = .N))]
  } else {
    f <- fold_override[run %chin% d$run]
    d[f, fold := i.fold, on = .(run, cb)]
    stopifnot(nrow(f) == nrow(d), !anyNA(d$fold))
  }
  depth <- depth_threshold(d, total_floor, total_cap, total_fraction)
  d <- merge(d, depth, by = "sample", sort = FALSE)
  d <- call_placements(rna_thresholds(d, primary_q, primary_floor, sparse), dominance, primary_beads)
  stopifnot(nrow(d) == sum(lengths(lapply(runs, `[[`, "cells"))),
    !anyDuplicated(d, by = c("run", "cb")), all(is.finite(d[accepted == TRUE, x])),
    all(is.finite(d[accepted == TRUE, y])), all(is.na(d[accepted == FALSE, x])),
    all(is.na(d[accepted == FALSE, y])))
  for (r in m$run[rows]) {
    z <- d[run == r]
    fwrite(z, file.path(out, "runs", paste0(r, "_audit.csv.gz")))
    fwrite(z[accepted == TRUE, .(sample, run, cb, x, y, primary_stable)],
      file.path(out, "runs", paste0(r, "_placements.csv")))
  }
  summaries[[s]] <- cbind(d[, .(eligible = .N, placed = sum(accepted),
    placement_rate = mean(accepted), rna_reference_unsupported = sum(rna_reference_status == "UNSUPPORTED"),
    beads_in_density_map = nrow(dm), beads_removed = sum(dm$density_removed),
    density_radius_um = median(m$eps_um[rows]))], depth)
}
fwrite(rbindlist(summaries), file.path(out, "sample_summary.csv"))
fwrite(m, file.path(out, "resolved_manifest.csv"))
fwrite(data.table(parameter = names(defaults), value = unname(defaults)), file.path(out, "parameters.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
writeLines("Complete", file.path(out, "DONE"))
