read_table <- function(path, ids = c("sample", "run", "cb", "sb", "group")) {
  if (!file.exists(path) || dir.exists(path)) stop("Missing file: ", path)
  columns <- names(fread(path, nrows = 0L, showProgress = FALSE))
  fread(path, colClasses = list(character = intersect(ids, columns)), showProgress = FALSE)
}

require_cols <- function(x, cols, label) {
  if (!all(cols %in% names(x))) stop(label, " missing columns: ", paste(setdiff(cols, names(x)), collapse = ", "))
}

valid_ids <- function(x) is.character(x) && !anyNA(x) && all(nzchar(x))

density_map <- function(beads, eps, cutoff = .9) {
  beads <- unique(beads[, .(sb, x, y)])
  if (anyDuplicated(beads$sb)) stop("Bead coordinates differ across runs sharing this puck")
  if (!nrow(beads)) stop("Cannot estimate density from an empty bead set")
  nn <- dbscan::frNN(as.matrix(beads[, .(x, y)]), eps, sort = FALSE)
  beads[, local_density := lengths(nn$id)]
  beads[, density_percentile := frank(local_density, ties.method = "average") / .N]
  beads[, density_removed := density_percentile > cutoff][]
}

cluster_one_cell <- function(links, eps) {
  bounds <- find_boundaries(links, eps)
  none <- list(bounds = bounds, nominal_n = 0L, candidates = NULL)
  mp <- bounds[["i2"]]
  if (mp < 1L) return(none)
  nominal <- cluster_table(links, db_clusters(links, eps, mp))
  none$nominal_n <- nrow(nominal)
  if (!nrow(nominal)) return(none)
  ref <- nominal[cluster_umi >= 5 & cluster_beads >= 2L & cluster_rms <= 100]
  if (!nrow(ref)) return(none)
  scale <- c(.8, .9, 1, 1.1, 1.2)
  grid <- rbind(
    data.table(axis = "scale", radius = eps * scale, minpts = mp * scale^2),
    data.table(axis = "evidence", radius = eps, minpts = mp * c(.75, .875, 1, 1.125, 1.25))
  )
  grid[, minpts := pmax(1L, as.integer(round(minpts)))]
  grid <- unique(grid)
  matches <- rbindlist(lapply(seq_len(nrow(grid)), function(j) {
    test <- cluster_table(links, db_clusters(links, grid$radius[j], grid$minpts[j]))
    match_clusters(ref, test)[, axis := grid$axis[j]][]
  }))
  metrics <- matches[, {
    k <- match(.BY$cluster_rank, ref$cluster_rank)
    ok <- matched & is.finite(test_x) & is.finite(test_y)
    shift <- sqrt((test_x[ok] - ref$cluster_x[k])^2 + (test_y[ok] - ref$cluster_y[k])^2)
    .(valid = mean(ok), jac = if (any(ok)) median(jaccard[ok]) else 0,
      shift = if (any(ok)) as.numeric(quantile(shift, .9)) else NA_real_)
  }, by = .(cluster_rank, axis)]
  metrics[, axis_pass := valid >= .5 & jac >= .8 & !is.na(shift) & shift <= 20]
  wide <- dcast(metrics, cluster_rank ~ axis, value.var = c("valid", "jac", "shift", "axis_pass"))
  ref <- merge(ref, wide, by = "cluster_rank", sort = FALSE)
  ref[, stable := axis_pass_evidence & axis_pass_scale]
  ref[, member_sbs := vapply(member_sbs, paste, collapse = ";", FUN.VALUE = character(1))]
  list(bounds = bounds, nominal_n = nrow(nominal), candidates = ref)
}

profile_run <- function(links, cells, density, eps, link_min, cores = 1L) {
  links <- links[, .(umi = sum(umi)), by = .(cb, sb, x, y)]
  raw <- links[, .(raw_spatial_umi = sum(umi)), by = cb]
  links <- links[umi >= link_min & !sb %chin% density[density_removed == TRUE, sb]]
  filtered <- links[, .(filtered_spatial_umi = sum(umi)), by = cb]
  by_cell <- split(links, by = "cb", keep.by = FALSE)
  run_one <- function(cb) cluster_one_cell(by_cell[[cb]], eps)
  profiles <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(cells, run_one, mc.cores = cores) else lapply(cells, run_one)
  if (any(vapply(profiles, inherits, logical(1), "try-error"))) stop("A clustering worker failed")
  candidates <- rbindlist(Map(function(z, id) {
    if (!is.null(z$candidates)) z$candidates[, cb := id][]
  }, profiles, cells))
  if (!nrow(candidates)) candidates <- data.table(cb = character(), cluster_rank = integer(),
    cluster_umi = numeric(), cluster_beads = integer(), cluster_rms = numeric(),
    cluster_x = numeric(), cluster_y = numeric(), stable = logical(), member_sbs = character())
  audit <- data.table(cb = cells, eps_um = eps,
    minPts1 = vapply(profiles, function(z) as.integer(z$bounds[["i1"]]), integer(1)),
    minPts2 = vapply(profiles, function(z) as.integer(z$bounds[["i2"]]), integer(1)),
    nominal_n = vapply(profiles, function(z) as.integer(z$nominal_n), integer(1)))
  audit[raw, raw_spatial_umi := i.raw_spatial_umi, on = "cb"]
  audit[filtered, filtered_spatial_umi := i.filtered_spatial_umi, on = "cb"]
  for (nm in c("raw_spatial_umi", "filtered_spatial_umi"))
    set(audit, which(is.na(audit[[nm]])), nm, 0)
  list(audit = audit, candidates = candidates)
}

summarize_candidates <- function(candidates) {
  c <- copy(candidates)
  setorder(c, run, cb, -cluster_umi, cluster_rank)
  if (!nrow(c)) return(data.table(run = character(), cb = character(), qualifying_n = integer(),
    primary_umi = numeric(), secondary_umi = numeric(), primary_beads = integer(),
    primary_x = numeric(), primary_y = numeric(), primary_rms = numeric(),
    primary_stable = logical(), stable_n = integer(), stable_ratio = numeric()))
  c[, {
    su <- cluster_umi[which(stable %in% TRUE)]
    .(qualifying_n = .N, primary_umi = cluster_umi[1L],
      secondary_umi = if (.N > 1L) as.numeric(cluster_umi[2L]) else NA_real_,
      primary_beads = cluster_beads[1L], primary_x = cluster_x[1L], primary_y = cluster_y[1L],
      primary_rms = cluster_rms[1L], primary_stable = stable[1L],
      stable_n = length(su), stable_ratio = if (length(su) > 1L) su[1L] / su[2L] else NA_real_)
  }, by = .(run, cb)]
}

depth_threshold <- function(d, floor = 30, cap = 100, fraction = .05, min_ref = 100L) {
  d <- copy(d)
  # The original provisional GOLD/accepted sets calibrate total UMI only.
  # Stability never directly rejects a final RNA-adjusted call.
  d[, gold := filtered_spatial_umi >= 10 & stable_n == 1L]
  d[, provisional := gold | (filtered_spatial_umi >= 10 & stable_n > 1L &
    !is.na(stable_ratio) & stable_ratio >= 3)]
  d[, {
    usable <- filtered_spatial_umi > 0
    ref <- if (sum(gold) >= min_ref) "provisional_GOLD" else if (sum(provisional) >= min_ref)
      "provisional_accepted_fallback" else if (any(usable)) "positive_UMI_fallback" else "no_positive_UMI"
    x <- switch(ref, provisional_GOLD = filtered_spatial_umi[gold],
      provisional_accepted_fallback = filtered_spatial_umi[provisional],
      positive_UMI_fallback = filtered_spatial_umi[usable], no_positive_UMI = numeric())
    .(depth_reference_source = ref, depth_reference = as.numeric(median(x)), depth_n = length(x),
      adaptive_threshold = if (length(x)) max(floor, min(cap, fraction * median(x))) else as.numeric(floor))
  }, by = sample]
}

rna_thresholds <- function(d, probability = .05, floor = 10, on_sparse = "fail") {
  d <- copy(d)
  d[, log_rna := log10(rna_umi)]
  d[, `:=`(rna_cutoff = NA_real_, rna_reference_n = 0L, rna_reference_status = "UNSUPPORTED")]
  for (s in unique(d$sample)) for (f in sort(unique(d[sample == s, fold]))) {
    train <- d[sample == s & fold != f & qualifying_n > 0, .(log_rna, primary_umi)]
    i <- which(d$sample == s & d$fold == f)
    b <- if (nrow(train) > 1L) unique(quantile(train$log_rna, seq(0, 1, .1))) else numeric()
    d[i, rna_reference_n := nrow(train)]
    if (length(b) > 1L) {
      train[, bin := cut(log_rna, b, include.lowest = TRUE)]
      curve <- train[, .(x = median(log_rna), y = as.numeric(quantile(log10(primary_umi), probability))),
        by = bin][order(x)]
      if (nrow(curve) >= 2L) {
        d[i, `:=`(rna_cutoff = ceiling(pmax(floor, 10^approx(curve$x, curve$y,
          xout = log_rna, rule = 2)$y)), rna_reference_status = "FIT")]
        next
      }
    }
    if (on_sparse == "floor") d[i, `:=`(rna_cutoff = ceiling(floor), rna_reference_status = "FLOOR_FALLBACK")]
  }
  d[, log_rna := NULL][]
}

call_placements <- function(d, dominance = 3, min_beads = 3L) {
  d <- copy(d)
  d[, dominance_ratio := fifelse(qualifying_n > 1L, primary_umi / pmax(secondary_umi, 1), NA_real_)]
  d[, status := fcase(
    qualifying_n == 0L, "NO_QUALIFYING_CLUSTER",
    filtered_spatial_umi < adaptive_threshold, "LOW_TOTAL_SPATIAL_UMI",
    qualifying_n > 1L & dominance_ratio < dominance, "AMBIGUOUS_COMPETITORS",
    primary_beads < min_beads, "TOO_FEW_PRIMARY_BEADS",
    is.na(rna_cutoff), "RNA_REFERENCE_UNSUPPORTED",
    primary_umi < rna_cutoff, "LOW_PRIMARY_UMI",
    default = "PLACED")]
  d[, accepted := status == "PLACED"]
  d[, `:=`(x = fifelse(accepted, primary_x, NA_real_), y = fifelse(accepted, primary_y, NA_real_))][]
}
