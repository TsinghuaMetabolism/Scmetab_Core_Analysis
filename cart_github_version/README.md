# CAR-T main-figure analysis

Computational methods for the CAR-T main figures (Figure 4b–g and Figure 5d–f).
This pre-publication snapshot contains source code and an empty parameter
template. Study matrices, annotations, fitted models, per-cell results and
dataset-specific parameter values are withheld until publication.

| Analysis | Code |
|---|---|
| Feature matching, original-TIC normalization and Bridge mean correction | `bounded_run_groups`, `tic_log1p`, `bridge_mean` in `CART_analysis.py` |
| Ordering-feature screen and rank-normal transform | `select_ordering` in `CART_analysis.py` |
| Ctr DDRTree reference, biological root and model export (4b–d) | `CART_trajectory.R` |
| Continuous feature dynamics, wild-bootstrap tests and joint heatmap selection (4e) | `fit_dynamics`, `joint_heatmap` in `CART_analysis.py` |
| Original-TIC spline curves with HC3 intervals (4f) | `write_spline_basis` in `CART_trajectory.R`; `fit_ols_hc3` in `CART_analysis.py` |
| Actual-day and same-run treatment comparisons (4g, 5d) | `abundance_statistics` in `CART_analysis.py` |
| Conditional treatment coordinates and frozen-tree projection (5e) | `reference_transform`, `solve_coordinates`, `freeze_tree`, `project_graph` in `CART_analysis.py` |
| Shared-axis pseudotime distribution and reflected KDE (5f) | `pseudotime_distribution` in `CART_analysis.py` |
| Main-panel rendering from external checkpoints | `run.py` and `plot_*.R` |

## Run the main panels

Install the Python dependencies with `python -m pip install -r requirements.txt`.
Rendering uses R with `ggplot2`, `data.table`, `jsonlite`, `ragg` and `mgcv`.
The computational R requirements are listed below. Use Python 3.11 or later.
Prepare the external bundle described in [DATA_ACCESS.md](DATA_ACCESS.md), then
fill the relevant fields in a local copy of `CART_parameters.example.json`.

```bash
python run.py render --input /path/to/bundle --output /path/to/new_run \
  --settings /path/to/CART_parameters.json --rscript /path/to/Rscript
```

Use `--panel figure4g` to run one panel. The output directory must be new and
outside the code directory and input bundle. It contains the figures, local
working tables, a run log and a manifest. The source bundle is read only.

This command recomputes Figure 4e selection, Figure 4f spline fits, Figure 4g/5d
abundance statistics and Figure 5f distributions. Figure 4b–d and Figure 5e
render supplied frozen trajectory/projection coordinates. These are checkpoint
operations; the library exposes the upstream algorithms separately. A complete
study run additionally requires the external feature matching, QC membership,
annotations, model settings and input bundle.

```bash
python run.py trajectory --input /path/to/ordering_input --output /path/to/new_reference \
  --settings /path/to/CART_parameters.json --rscript /path/to/Rscript
python run.py basis --input /path/to/spline_input --output /path/to/new_basis \
  --settings /path/to/CART_parameters.json --rscript /path/to/Rscript
```

`trajectory` fits and exports the Ctr reference; `basis` exports the spline
design. Their input contracts are described below. Projection functions in
`CART_analysis.py` accept the frozen transform, loadings, graph and query arrays.

## Inputs and settings

Python functions accept arrays and data frames; their docstrings specify the
required columns. Matrices use **cells × features**, with the same cell order as
the metadata. Original non-detections are zero; features not measured in a run
are `NaN`. Each run retains its original TIC denominator. Cross-run feature
identities must be aligned before Bridge correction. A feature group contains
at most one feature from each run and satisfies the supplied ppm span bound.

Copy `CART_parameters.example.json` to a local settings file and supply the
study-specific values. The template contains `null` placeholders and is not an
executable study configuration. Keep the input data and completed settings
outside the repository.

For a single Ctr run, normalize the raw intensities directly. For the combined
Ctr/Myr analysis, align the feature axes and apply the corresponding Bridge
baseline mean shifts before the reference and projection steps. Ctr reference
ordering uses finite features; normalized intensity remains the response for
abundance/dynamic analysis, while the rank-normal matrix defines geometry.

The R fitting command reads `ordering.tsv.gz` (**features × cells**, first
column feature ID) and `cells.tsv` (first column cell ID, numeric `round` column)
from a local input directory:

```bash
Rscript CART_trajectory.R fit INPUT_DIR OUTPUT_DIR SETTINGS.json
```

It writes the rooted reference, coordinates, learned nodes/edges, loadings and
native tree projections. Spline construction reads `pseudotime.tsv` and
`prediction_grid.tsv`, each with a `t` column:

```bash
Rscript CART_trajectory.R basis INPUT_DIR OUTPUT_DIR SETTINGS.json
```

Both exported design matrices include an intercept. Fit the dynamic response
with `fit_dynamics`. For the five target curves, use `fit_ols_hc3` on the same
saved basis with original-TIC intensities; calculate pointwise mean-curve
intervals from the returned covariance.

Treatment projection reuses the reference feature transform, loadings, nodes,
root and reference-supported leaf extensions. Compute Ctr conditional
coordinates using the same feature panel available in Myr before comparing
their positions. Retain native Ctr coordinates for the Figure 5e background and
use conditional graph positions for the Figure 5f comparison.

## Statistical conventions

Dynamic tests use a null-imposed Mammen wild bootstrap and BH adjustment over
all finite, variable features in the corresponding reference analysis.
Cross-donor heatmap rows require concordant fitted endpoint directions and the
specified number of donors passing the dynamic-q threshold. Mean pairwise
correlation determines row order. Prespecified SM rows receive display priority
after satisfying the same inclusion rule.

Abundance tests use two-sided asymptotic Mann–Whitney U with tie and continuity
correction. Each panel forms one BH family across all prespecified donors,
targets and group pairs. Group means and SEM include zero observations.
Pseudotime distributions use one shared frozen axis; each group is not rescaled
separately. These analyses describe cell-level distributions and associations
conditional on the reference trajectory.

## Software

Python: NumPy, pandas and SciPy (NumPy 2 or later is required for `trapezoid`).
Reading original H5AD inputs separately requires AnnData. R: Monocle 2, DDRTree, VGAM, Biobase,
igraph, mgcv and jsonlite. The study environment used Monocle 2.34.0, DDRTree
0.1.5 and igraph 2.0.3. No notebook execution output or study results are
embedded in these files.

## Validation

The release functions were checked locally against the accepted main-panel
results for abundance tests, heatmap selection/order, reflected KDE, spline
bases and sampled treatment projections. A Ctr reference refit matched the
saved reference within floating-point tolerance. All nine panel entry points
were exercised on the external study bundle. Study-specific validation outputs
and figures remain outside this code snapshot.
