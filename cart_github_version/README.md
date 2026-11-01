# CAR-T metabolic trajectories

Computational analysis for CAR-T Figure 4b–g and Figure 5d–f. The package follows the NSCLC release layout: a root launcher, documentation, workflow orchestration, individual panel modules and upstream scientific methods.

## Repository layout

```text
cart_github_version/
  README.md
  run.py
  requirements.txt
  environment.json
  config/
    CART_parameters.example.json
  docs/
    DATA_ACCESS.md
    METHODS.md
    code_manifest.json
  scripts/
    run_workflow.py
  workflow/
    provenance/panels.json
    scripts/
      panel_statistics.py
      plot_common.R
    upstream/
      preprocessing/preprocess.py
      trajectory/fit_reference.R
      dynamics/fit_dynamics.py
      projection/project_reference.py
    panel_4b-T1_trajectory/scripts/
    panel_4c-T3_trajectory/scripts/
    panel_4d-Batch03_trajectory/scripts/
    panel_4e-shared_metabolic_dynamics/scripts/
    panel_4f-SM_metabolic_curves/scripts/
    panel_4g-SM_day_abundance/scripts/
    panel_5d-Bridge_SM_comparison/scripts/
    panel_5e-Myr_projection/scripts/
    panel_5f-pseudotime_distribution/scripts/
```

Each panel directory contains its R plotting implementation (`plot.R`) and a Python entry point (`plot.py`). Shared algorithms and rendering helpers have one maintained implementation.

## Analysis modules

| Module | Code | Purpose |
|---|---|---|
| Preprocessing | `workflow/upstream/preprocessing/` | Feature matching, original-TIC normalization, Bridge mean correction and ordering-feature selection |
| Ctr reference | `workflow/upstream/trajectory/` | DDRTree fitting, biological rooting, model export and cubic-regression spline basis |
| Dynamic metabolites | `workflow/upstream/dynamics/` | Continuous spline fits, HC3 intervals, Mammen bootstrap tests and full-family BH |
| Fixed-reference projection | `workflow/upstream/projection/` | Reference transform, conditional coordinates and frozen-tree projection |
| Main-panel statistics | `workflow/scripts/panel_statistics.py` | Abundance comparisons, shared-donor heatmap selection and pseudotime distributions |
| Figure reconstruction | `workflow/panel_*/scripts/` | Nine main-figure panels from the external study checkpoints |

## Installation

Use Python 3.11 or later and install dependencies with `python -m pip install -r requirements.txt`. Rendering requires R with `ggplot2`, `data.table`, `jsonlite`, `ragg` and `mgcv`. Trajectory fitting also requires Monocle 2, DDRTree, VGAM, Biobase and compatible igraph. Tested environments are recorded in [environment.json](environment.json).

## Data and settings

Study matrices, annotations, fitted models, per-cell results and completed parameter values are withheld until publication. The exact external bundle is specified in [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md). Copy `config/CART_parameters.example.json` outside the repository and fill the required null fields with the approved study configuration.

## Run

Render all main panels:

```bash
python run.py render --input /path/to/bundle --output /path/to/new_run \
  --settings /path/to/CART_parameters.json --rscript /path/to/Rscript
```

Render a single panel using the root launcher or its own entry point:

```bash
python run.py render --panel figure4g --input /path/to/bundle \
  --output /path/to/new_run --settings /path/to/CART_parameters.json
python workflow/panel_4g-SM_day_abundance/scripts/plot.py \
  --input /path/to/bundle --output /path/to/new_run \
  --settings /path/to/CART_parameters.json
```

The output directory must be new and outside the code directory and input bundle. It contains figures, local working tables, a run log and a manifest. Source data are read only. Figure 4e selection, Figure 4f spline fits, Figure 4g/5d abundance statistics and Figure 5f distributions are recomputed; Figure 4b–d and Figure 5e render supplied frozen coordinates.

Fit a Ctr reference or export a spline design:

```bash
python run.py trajectory --input /path/to/ordering_input --output /path/to/new_reference \
  --settings /path/to/CART_parameters.json --rscript /path/to/Rscript
python run.py basis --input /path/to/spline_input --output /path/to/new_basis \
  --settings /path/to/CART_parameters.json --rscript /path/to/Rscript
```

The upstream Python modules expose preprocessing, dynamic fitting and projection functions for the external study inputs. These are checkpoint operations; a complete study run additionally requires the approved QC membership, feature mapping, annotations and model configuration. See [docs/METHODS.md](docs/METHODS.md) for model and statistical conventions.

## Validation

All nine panel entry points were exercised locally. Abundance tests, heatmap selection/order, fitted curves and pseudotime distributions matched the accepted results. Reference fitting and sampled treatment projections were checked numerically. Study data and validation outputs remain outside this source package.
