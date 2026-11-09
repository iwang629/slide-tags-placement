# Weighted DBSCAN boundary, centroid and stability-matching primitives.
# Preserves the nominal clustering used by the analyzed coordinate release.
db_clusters <- function(dl, eps, minpts) {
  if (is.null(dl) || !nrow(dl) || !is.finite(minpts) || minpts < 1L) return(integer())
  dbscan::dbscan(dl[, .(x, y)], eps = eps, minPts = as.integer(minpts),
    weights = as.numeric(dl$umi), borderPoints = TRUE)$cluster
}

iestimate <- function(dl, eps) {
  if (is.null(dl) || !nrow(dl)) return(0L)
  w <- as.numeric(dl$umi)
  qry <- as.matrix(dl[w == max(w), .(x, y)])
  nn <- dbscan::frNN(as.matrix(dl[, .(x, y)]), eps, query = qry, sort = FALSE)
  as.integer(floor(max(vapply(nn$id, function(ii) sum(w[ii]), numeric(1)))))
}

ipercluster <- function(dl, v, eps) {
  ids <- sort(unique(v[v > 0L]))
  vapply(ids, function(cid) {
    s <- dl[v == cid]; w <- as.numeric(s$umi)
    nn <- dbscan::frNN(as.matrix(s[, .(x, y)]), eps, sort = FALSE)
    vals <- vapply(seq_along(nn$id), function(i) sum(w[nn$id[[i]]]) + w[i], numeric(1))
    as.integer(floor(max(vals)))
  }, integer(1))
}

find_boundaries <- function(dl, eps) {
  if (is.null(dl) || !nrow(dl)) return(c(i1 = 0L, i2 = 0L))
  if (nrow(dl) == 1L) return(c(i1 = as.integer(sum(dl$umi)), i2 = 0L))
  i1 <- tryCatch({
    i <- iestimate(dl, eps); v <- db_clusters(dl, eps, i)
    stopifnot(length(v) > 0L, max(v) >= 1L)
    i <- max(ipercluster(dl, v, eps)); v <- db_clusters(dl, eps, i)
    stopifnot(length(v) > 0L, max(v) >= 1L)
    while (max(v) > 0L) {i <- i + 1L; v <- db_clusters(dl, eps, i)}
    as.integer(i - 1L)
  }, error = function(e) 0L)
  i2 <- tryCatch((function() {
    i <- iestimate(dl, eps); v <- db_clusters(dl, eps, i)
    stopifnot(length(v) > 0L, max(v) >= 1L)
    seen <- integer()
    repeat {
      seen <- c(seen, i)
      if (any(v != 1L)) {
        i <- iestimate(dl[v != 1L], eps); v <- db_clusters(dl, eps, i)
      } else return(0L)
      if (length(v) && max(v) >= 2L) {
        ip <- sort(ipercluster(dl, v, eps)); if (length(ip) < 2L) return(0L)
        i <- ip[length(ip) - 1L]; v <- db_clusters(dl, eps, i)
      }
      if ((length(v) && max(v) >= 2L) || i %in% seen) break
    }
    if (!length(v)) return(0L)
    if (max(v) < 2L) {
      while (max(v) < 2L && i >= 1L) {i <- i - 1L; v <- db_clusters(dl, eps, i)}
      return(as.integer(i))
    }
    while (max(v) >= 2L) {i <- i + 1L; v <- db_clusters(dl, eps, i)}
    as.integer(i - 1L)
  })(), error = function(e) 0L)
  c(i1 = i1, i2 = as.integer(min(i2, i1)))
}

cluster_table <- function(dl, cl) {
  ids <- sort(unique(cl[cl > 0L])); if (!length(ids)) return(data.table())
  z <- rbindlist(lapply(ids, function(cid) {
    d <- dl[cl == cid]; cx <- weighted.mean(d$x, d$umi); cy <- weighted.mean(d$y, d$umi)
    radial <- sqrt((d$x - cx)^2 + (d$y - cy)^2)
    data.table(cluster_id = cid, cluster_umi = sum(d$umi), cluster_beads = nrow(d),
      cluster_x = cx, cluster_y = cy,
      cluster_rms = sqrt(weighted.mean(radial^2, d$umi)),
      member_sbs = list(as.character(d$sb)))
  }))
  setorder(z, -cluster_umi, cluster_id); z[, cluster_rank := seq_len(.N)]; z[]
}

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (!u) NA_real_ else length(intersect(a, b)) / u
}

match_clusters <- function(ref, test) {
  if (!nrow(ref)) return(data.table())
  if (!nrow(test)) return(ref[, .(cluster_rank, matched = FALSE, jaccard = 0,
    test_x = NA_real_, test_y = NA_real_)])
  pairs <- CJ(ri = seq_len(nrow(ref)), ti = seq_len(nrow(test)))
  pairs[, jac := mapply(function(i, j) jaccard(ref$member_sbs[[i]],
    test$member_sbs[[j]]), ri, ti)]
  setorder(pairs, -jac, ri, ti)
  used_r <- integer(); used_t <- integer(); chosen <- list()
  for (k in seq_len(nrow(pairs))) {
    p <- pairs[k]
    if (p$ri %in% used_r || p$ti %in% used_t || !is.finite(p$jac) || p$jac <= 0) next
    chosen[[length(chosen) + 1L]] <- p
    used_r <- c(used_r, p$ri); used_t <- c(used_t, p$ti)
  }
  z <- ref[, .(cluster_rank)]
  z[, `:=`(matched = FALSE, jaccard = 0, test_x = NA_real_, test_y = NA_real_)]
  if (length(chosen)) {
    ch <- rbindlist(chosen)
    for (k in seq_len(nrow(ch))) {
      i <- ch$ri[k]; j <- ch$ti[k]
      z[cluster_rank == ref$cluster_rank[i], `:=`(matched = TRUE,
        jaccard = ch$jac[k], test_x = test$cluster_x[j], test_y = test$cluster_y[j])]
    }
  }
  z[]
}
