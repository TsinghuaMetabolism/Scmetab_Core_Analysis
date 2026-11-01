# NSCLC metabolic trajectories

Computational analysis for NSCLC Figure 3, panels b–i. Clinical information and experimental viability data are outside this repository.

## Analysis modules

| Module | Code | Purpose |
|---|---|---|
| P59 preprocessing | `workflow/upstream/P59/02_analysis/P59_final_analysis_20260919/scripts/` | Recorded-TIC normalization, dual-anchor correction, saved-ordering DDRTree refit |
| P39 | `workflow/upstream/P39/02_analysis/` | Endpoint-consistent ordering selection and current downstream Gaussian GAM |
| P41 | `workflow/upstream/P41/02_analysis/` | Import, feature matching, integration, reference trajectory, projection and GAM |
| P59 dynamic metabolites | `workflow/upstream/P59/02_analysis/P59_choline_25_GAM_full1928_20261005/` | Full-feature GAM with and without D0 |
| Shared dynamics / pathways | `workflow/panel_f-common_pathway_enrichment/` | Shared-feature selection, KEGG ORA and independent validation |
| Figure reconstruction | `workflow/scripts/rebuild.py` | English HTML and panels b–i |
| Distribution | `workflow/panel_i-P59_pseudotime_distribution/scripts/` | Fixed-reference distributions and reflected KDE |

## Installation

Python 3.12 and R are required. Install Python dependencies with `python -m pip install -r requirements-upstream.txt`. Tested versions are in `environment.json`. R trajectory analysis requires Monocle2, DDRTree, VGAM, Biobase, data.table and compatible igraph; GAM requires mgcv. Use a compatible existing environment rather than blindly updating Monocle dependencies.

## Data availability

Expression/intensity matrices are not distributed in this repository. Request the approved research data bundle from the study authors. The data-access contact and terms must be added by the authors before public release. See [data requirements](docs/DATA_ACCESS.md).

## Validated execution

The separately supplied research bundle must retain its original relative layout and `provenance/source_manifest.csv`. Its path is passed explicitly; no author-specific computer path is required.

```bash
python run.py audit --package /path/to/approved_bundle --output /path/to/new_run
python run.py redraw --package /path/to/approved_bundle --output /path/to/new_run
python run.py p59-raw --package /path/to/approved_bundle --output /path/to/new_run
python run.py p59-trajectory --package /path/to/approved_bundle --output /path/to/new_run --rscript /path/to/Rscript
python run.py p59-gam --package /path/to/approved_bundle --output /path/to/new_run --rscript /path/to/Rscript
```

Outputs are written into a separate run directory, including `package/index.html`, figures, numerical tables and JSON reports. The source bundle is not modified. The validated launcher currently executes the approved bundle's scientific sources; the `workflow/` tree exposes those sources for review. Do not treat editing `workflow/` as changing an already supplied bundle automatically.

## Reproduction scope

Figure replay, P59 raw-to-corrected recomputation, saved-ordering trajectory refit and current-projection GAM refit were independently tested. These are individual checkpoint-level operations, not a connected fresh raw-to-all-panels pipeline. P59 raw recomputation retains frozen QC and feature matching. P41 original B2/B3 H5AD and original annotation CSV remain missing from the available archive. Current query projections and shared-donor curves remain supplied checkpoints.

A matrix alone cannot reproduce the published panels: cell metadata, feature annotation, reference definitions and the documented checkpoints/resources are also required. Archived P39/P41 scripts are scientific source modules and may still contain historical path assumptions; use the validated commands above for the supported reproduction levels. See [methods and release boundaries](docs/METHODS.md).

Repository organization follows the concise task-oriented style of https://github.com/guomics-lab/PTV-1; no code was copied from that repository. No license is assigned here without the authors' decision.
