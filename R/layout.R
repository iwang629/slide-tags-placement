# Optional downstream utilities. They never infer tissue groups or change
# placement decisions. Positions are in micrometers.
fit_rigid <- function(x, target, min_anchors = 20L) {
  if (nrow(x) < min_anchors || nrow(target) != nrow(x) ||
      !all(is.finite(x)) || !all(is.finite(target))) stop("Insufficient valid alignment anchors")
  fit <- function(keep) {
    a <- x[keep, , drop = FALSE]; b <- target[keep, , drop = FALSE]
    xc <- colMeans(a); yc <- colMeans(b)
    cross <- crossprod(sweep(a, 2, xc), sweep(b, 2, yc))
    sv <- svd(cross)
    if (sv$d[1L] <= 0 || sv$d[2L] <= 1e-10 * sv$d[1L]) stop("Degenerate alignment anchors")
    q <- sv$u %*% t(sv$v)
    shift <- yc - as.vector(xc %*% q)
    err <- sqrt(rowSums((sweep(x %*% q, 2, shift, "+") - target)^2))
    list(q = q, shift = shift, error = err)
  }
  keep <- rep(TRUE, nrow(x))
  for (i in seq_len(15L)) {
    err <- fit(keep)$error
    next_keep <- err <= max(5, median(err) + 6 * mad(err))
    if (sum(next_keep) < min_anchors) stop("Too few robust alignment inliers")
    if (identical(keep, next_keep)) break
    keep <- next_keep
  }
  c(fit(keep), list(keep = keep))
}

organize_placements <- function(placements, groups, anchors = NULL) {
  p <- copy(as.data.table(placements))
  g <- as.data.table(groups)
  keys <- c("run", "cb")
  stopifnot(all(c("sample", keys, "x", "y") %in% names(p)),
    all(c(keys, "group") %in% names(g)), !anyDuplicated(p, by = keys),
    !anyDuplicated(g, by = keys), all(is.finite(p$x)), all(is.finite(p$y)))
  p[, `:=`(group = NA_character_, layout_x = NA_real_, layout_y = NA_real_)]
  p[g, group := as.character(i.group), on = keys]
  p[!is.na(group) & !nzchar(group), group := NA_character_]
  empty <- list(coordinates = p, transforms = data.table())
  if (is.null(anchors) || !any(!is.na(p$group))) return(empty)
  a <- as.data.table(anchors)
  stopifnot(all(c(keys, "group", "target_x", "target_y") %in% names(a)),
    !anyDuplicated(a, by = keys), all(is.finite(a$target_x)), all(is.finite(a$target_y)))
  ref <- merge(p[!is.na(group)], a, by = c(keys, "group"))
  groups_present <- unique(p[!is.na(group), .(sample, group)])
  transforms <- rbindlist(lapply(seq_len(nrow(groups_present)), function(i) {
    ss <- groups_present$sample[i]; gg <- groups_present$group[i]
    rr <- ref[sample == ss & group == gg]
    f <- fit_rigid(as.matrix(rr[, .(x, y)]), as.matrix(rr[, .(target_x, target_y)]))
    ix <- which(p$sample == ss & p$group == gg)
    xy <- sweep(as.matrix(p[ix, .(x, y)]) %*% f$q, 2, f$shift, "+")
    set(p, ix, "layout_x", xy[, 1L]); set(p, ix, "layout_y", xy[, 2L])
    data.table(sample = ss, group = gg, anchors = nrow(rr), inliers = sum(f$keep),
      rmse_um = sqrt(mean(f$error[f$keep]^2)), reflection = det(f$q) < 0,
      q11 = f$q[1, 1], q12 = f$q[1, 2], q21 = f$q[2, 1], q22 = f$q[2, 2],
      shift_x = f$shift[1], shift_y = f$shift[2])
  }))
  list(coordinates = p, transforms = transforms)
}

attach_seurat <- function(object, coordinates, run, barcode_map = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject")
  cells <- colnames(object)
  p <- as.data.table(coordinates)
  stopifnot(all(c("run", "cb", "x", "y") %in% names(p)), !anyDuplicated(p, by = c("run", "cb")))
  p <- p[which(p$run == run)]
  if (!all(is.finite(p$x) & is.finite(p$y))) stop("Use the accepted placements table only")
  if (!is.null(barcode_map)) {
    stopifnot(all(c("cb", "seurat_cell") %in% names(barcode_map)),
      !anyDuplicated(barcode_map$cb), !anyDuplicated(barcode_map$seurat_cell))
    p <- p[which(p$cb %in% barcode_map$cb)]
    p[, seurat_cell := barcode_map$seurat_cell[match(cb, barcode_map$cb)]]
  } else p[, seurat_cell := cb]
  p <- p[seurat_cell %in% cells]
  for (nm in setdiff(c("group", "layout_x", "layout_y"), names(p))) p[, (nm) := NA]
  i <- match(cells, p$seurat_cell)
  meta <- data.frame(spatial_group = as.character(p$group[i]), spatial_puck_x = p$x[i],
    spatial_puck_y = p$y[i], spatial_layout_x = as.numeric(p$layout_x[i]),
    spatial_layout_y = as.numeric(p$layout_y[i]), spatial_placed = !is.na(i),
    spatial_grouped = !is.na(p$group[i]), row.names = cells)
  names_to_add <- c(names(meta), "spatial_puck", "spatial_layout")
  if (any(names_to_add %in% c(names(object), colnames(object[[]]))))
    stop("Spatial fields already exist; attach to a clean copy")
  object <- SeuratObject::AddMetaData(object, meta)
  add <- function(obj, nm, x, y, prefix) {
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) return(obj)
    mat <- cbind(x[ok], y[ok])
    rownames(mat) <- cells[ok]
    colnames(mat) <- paste0(prefix, 1:2)
    obj[[nm]] <- SeuratObject::CreateDimReducObject(embeddings = mat,
      key = prefix, assay = SeuratObject::DefaultAssay(obj))
    obj
  }
  object <- add(object, "spatial_puck", meta$spatial_puck_x, meta$spatial_puck_y, "PUCK_")
  add(object, "spatial_layout", meta$spatial_layout_x, meta$spatial_layout_y, "LAYOUT_")
}
