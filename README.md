# RNA-adjusted Slide-tags spatial placement

Assign a puck coordinate to each called nucleus from its cell–bead links and
raw RNA UMI count. This implementation uses the strongest qualifying nominal
cluster and a within-sample RNA-adjusted support threshold. **Stability is a
QC flag, not a placement filter**. No legacy Coord1/Coord2 file, external training
cohort, tissue labels, or technical replicate is required.

## Try it with synthetic data

Requires R ≥4.3 and the R packages `data.table`, `dbscan`, and `R.utils`.
From the repository root:

```sh
Rscript -e 'install.packages(c("data.table", "dbscan", "R.utils"), repos="https://cloud.r-project.org")'
Rscript examples/make_demo.R demo_inputs
Rscript run_placement.R demo_inputs/manifest.csv demo_results
Rscript tests/smoke.R
```

The example contains two synthetic lanes on one puck. The demo and output
directories must not already exist. See `demo_results/sample_summary.csv` and
`demo_results/runs/lane1_placements.csv`.

## Use your data

Create a CSV manifest with one row per lane. `sample` is a **physical puck
group**: reuse it only when lanes share the same bead coordinates. Give each
run a unique ID; the cell key is **(run, cb)**, not the barcode alone. Paths are
absolute or relative to the manifest's folder. Use identifiers made of letters,
digits, periods, underscores and hyphens, beginning with a letter or digit.

```csv
sample,run,links,cells,rna,eps_um
puckA,lane1,lane1_links.csv.gz,lane1_cells.txt,lane1_rna.csv,49.5
puckA,lane2,lane2_links.csv.gz,lane2_cells.txt,lane2_rna.csv,49.5
```

| Input | Required contents |
| --- | --- |
| `links` | CSV or CSV.gz: `cb,sb,umi,x,y`. UMI is a positive integer count per cell–bead link. If a link is split across rows, counts are added. Bead coordinates are physical **µm** and must agree across lanes. |
| `cells` | One eligible barcode per line, including cells without usable links. Defines the placement-rate denominator. |
| `rna` | CSV: `cb,rna_umi`. Positive integer **raw RNA count sum** per eligible cell, not normalized expression, gene count, spatial UMI, or read depth. |
| `eps_um` | Nominal DBSCAN radius in µm for that lane. If omitted or blank, provide `eps_metadata` pointing to a one-row CSV with an `eps_um` column. |

Use exactly the same barcode strings across inputs. For a Seurat RNA assay,
sum the raw counts or use `nCount_RNA` if it refers to those counts. Normalize
barcode suffixes explicitly before running if your tools use different
conventions, and check that the mapping is one-to-one.

```sh
Rscript run_placement.R manifest.csv results --cores 4
Rscript run_placement.R --help
```

On Linux/macOS, cells can be processed in parallel; on Windows the pipeline
runs serially. It holds one puck's lanes in memory at a time. Treat
`DONE` as the completion marker and use a fresh directory for every rerun.

### Start from `SBcounts.h5` (optional)

Use a [Macosko-Pipelines Slide-tags](https://github.com/MacoskoLab/Macosko-Pipelines/tree/main/slide-tags)
HDF5 file that contains cell/bead observations **and** the embedded
bead-coordinate table, plus your RNA called-cell whitelist. A puck-coordinate
file alone is not enough. The HDF5 adapter requires `RANN` (CRAN) and
`rhdf5` (Bioconductor):

```sh
Rscript -e 'install.packages(c("RANN", "BiocManager"), repos="https://cloud.r-project.org"); BiocManager::install("rhdf5", ask=FALSE)'
Rscript prepare_sbcounts.R SBcounts.h5 called_cells.txt prepared_lane
```

Use `prepared_lane/links.csv.gz`, `prepared_lane/cells.txt`, and
`prepared_lane/metadata.csv` for `links`, `cells`, and `eps_metadata`,
respectively, in the manifest. Supply RNA counts separately. This adapter strips
10x numeric suffixes from the whitelist; make the RNA barcodes match. It
counts distinct spatial UMI identities after resolving competing beads by a
unique maximum read count. The default radius is the median distance to
the 36th returned neighbor (self included), using up to 10,000 seeded bead
queries. The optional fourth argument changes this index. This is an input
calibration choice, not a universal physical constant.

The legacy reconstruction route is also supported: supply `matrix.csv.gz`
as `links` (extra columns are ignored), and supply the intended cell list
and saved radius separately. For historical comparison, match the original
cell cohort, original radius and fold assignments. A newly estimated radius
or a broader whitelist can change coordinates, thresholds, and the denominator.

## Method and defaults

1. **Bead evidence:** Collapse duplicate cell–bead rows and retain links with
   ≥2 spatial UMIs. The density mask is built **before** the link filter from
   distinct observed beads in all lanes on that puck. Count neighboring beads
   within the median lane radius, rank counts with average ties, and remove
   beads above the 90th percentile of ranks. Ties can make the removed fraction
   differ from 10%. Do not compute density across different pucks.
2. **Nominal clusters:** Use the supplied run radius and UMI-weighted DBSCAN.
   A bead contributes its retained UMIs to neighborhood support. The
   per-cell `minPts2` is the searched upper boundary where a second cluster
   can still be detected. If the boundary is not found, the implementation does
   not substitute a different `minPts`. Qualifying nominal clusters have
   ≥5 UMI, ≥2 beads and weighted RMS radius ≤100 µm. The centroid is
   `sum(UMI × coordinate) / sum(UMI)`; weighted RMS is the square root of the
   weighted mean squared bead distance from that centroid.
3. **Stability and total evidence:** Re-run clustering on two sensitivity axes:
   radius multipliers 0.8–1.2 with `minPts` multiplied by the squared radius
   multiplier, and fixed radius with `minPts` multipliers 0.75–1.25. On both axes,
   the stability flag requires ≥50% successful matches, median bead-set
   Jaccard ≥0.8 and 90th-percentile centroid shift ≤20 µm. The historical
   *provisional* stable calls calibrate the puck's total filtered spatial UMI
   floor: `max(30, min(100, 0.05 × reference_median))`. Use ≥100 provisional
   GOLD cells (one stable candidate), else ≥100 provisional accepted cells
   (including cells with a ≥3-fold dominant stable candidate), else all
   positive-total cells. Provisional calls require ≥10 total UMI; when there are
   no positive totals, use the floor. Stability is **not** a direct rejection
   rule for the final placement.
4. **Strongest-first selection:** Rank the basic-qualifying nominal clusters
   by their spatial UMIs, with cluster rank breaking ties. If there is a
   second candidate, require the strongest to have ≥3× its UMI, even if that
   second candidate is unstable. Never fall back to a weaker cluster.
5. **RNA adjustment:** In each puck and held-out fold, use the other four folds'
   cells with a basic-qualifying primary cluster. Bin log RNA UMI by decile,
   calculate the 5th percentile of log primary-cluster spatial UMI in each bin,
   interpolate at the held-out cell's RNA UMI, hold endpoints constant, and
   round up with a floor of 10 cluster UMI. Require ≥3 retained beads and a
   primary-cluster UMI at or above this cutoff, in addition to the **separate**
   total spatial UMI cutoff. A passing cell is placed at the nominal primary
   centroid. There is no hard K, score, or SNR gate.

Seeded folds are stable for an unchanged set of run/barcode keys; use
`--folds folds.csv` with columns `run,cb,fold` (values 1–5) to replay an
archived assignment. Default behavior for a fold that cannot fit at least two
distinct RNA bins is to mark it unsupported. Explicitly pass
`--sparse-rna-reference floor` for a labeled floor-only fallback. This
fallback is for adaptation, not part of the historical fitted analysis.

## Outputs and interpretation

- `runs/*_placements.csv`: accepted coordinates only (`sample,run,cb,x,y`
  and the primary stability flag).
- `runs/*_audit.csv.gz`: every eligible cell, candidate support, `minPts`
  boundaries, fold, thresholds, status and first failure reason. Failed
  final coordinates are missing; candidate coordinates remain in the audit.
- `runs/*_candidates.csv.gz`, `density/*.csv.gz`: qualifying candidates and
  bead density masks. `sample_summary.csv`: eligible/placed counts and
  depth-reference source. `parameters.csv` and `session_info.txt` record settings
  and software versions.

For new pucks, inspect the density mask and failed/placed positions against
known tissue landmarks. A placement fraction, high bead support, or agreement
with another caller does not establish spatial accuracy. The RNA adjustment
addresses a difference in evidence depth; it does not correct all cell-type
selection bias. Folds calibrate a cutoff, not a probability of correct placement.

## Optional manual layout and Seurat attachment

See [manual grouping and Seurat](docs/SEURAT.md). Tissue grouping is separate
from coordinate calling; the placement script does not infer anatomy or
assign unlabeled cells to the nearest group.

This repository contains the placement implementation and synthetic examples
only. The upstream reconstruction workflow is separate:
[Macosko-Pipelines reconstruction](https://github.com/MacoskoLab/Macosko-Pipelines/tree/main/reconstruction).
See [LICENSE](LICENSE) for this package.
