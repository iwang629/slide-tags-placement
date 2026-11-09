# Manual tissue groups and Seurat

After reviewing raw puck coordinates, make an optional `groups.csv` keyed by
`run,cb` with a `group` column. Leave unassigned cells blank or NA. The same
group name can appear on different pucks. Use distinct labels for distinct
physical pieces on one puck.

If matching a previously reviewed layout, also create `anchors.csv` with
`run,cb,group,target_x,target_y`. The target is the cell's coordinate in the
previous layout, in µm. Use at least 20 reliable, non-collinear shared anchors
per piece. For example, retain cells whose manual group is unchanged and whose
new raw puck coordinate is within 50 µm of the prior raw coordinate.

```r
library(data.table)
source("R/placement.R")
source("R/layout.R")

p <- rbindlist(lapply(Sys.glob("results/runs/*_placements.csv"), read_table))
groups <- read_table("groups.csv")
anchors <- read_table("anchors.csv")
layout <- organize_placements(p, groups, anchors)
fwrite(layout$coordinates, "reviewed_coordinates.csv")
fwrite(layout$transforms, "transforms.csv")
```

The transform uses iterative trimming and fits rotation or reflection plus
translation. There is no scaling, per-cell displacement correction, or
nearest-neighbor group assignment. Check the returned transform diagnostics.
With no previous layout, use `organize_placements(p, groups)` to keep manual
groups without generating organized coordinates.

To attach only matching coordinates to a new copy of one run's Seurat object:

```r
obj <- readRDS("lane1_original.rds")
# Omit barcode_map if the IDs match exactly. Otherwise, supply a unique CSV
# mapping columns cb,seurat_cell (e.g. AACG...,AACG...-1).
mapped <- attach_seurat(obj, layout$coordinates, "lane1")
stopifnot(identical(colnames(obj), colnames(mapped)))
saveRDS(mapped, "lane1_with_spatial.rds")
```

Raw positions are in the `spatial_puck` reduction; organized coordinates, if
available, are in `spatial_layout`. RNA cells without placements are retained
with missing spatial coordinates. Placements absent from the target object
are not attached. The reductions are coordinate embeddings, not tissue images.
