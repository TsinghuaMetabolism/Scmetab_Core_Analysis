# HSC Figure 2 analysis code

This repository contains the analysis code and input contract for the HSC single-cell metabolomics Figure 2 study. It contains no expression matrices, cell-level metadata, fitted results, or instrument files.

## Scope

Given the requested **processed TIC-normalized log1p intensity matrix** and matching metadata, `run.py fit` fits HSC path-specific and MPP subset Tobit models, exports statistics and curves, then builds an offline English report for panels b–l and n–o. Supplied trajectory and UMAP coordinates are plotted; they are not estimated here.

The code does not process instrument RAW files, integrate peaks across batches, perform cell QC, estimate a new trajectory, or infer metabolite structure from a mass match. A bare raw intensity matrix is **not** an accepted model input. `code/prepare_input.py` can perform the stated TIC scaling and log1p transform on an *already integrated* nonnegative intensity matrix if matching per-cell TICs and the prespecified prefilter TIC reference are supplied. It does not perform the other upstream stages.

## Obtain the data

The expression matrices and cell-level data are available upon request from the study authors. See [DATA_ACCESS.md](DATA_ACCESS.md) for the exact external bundle layout and schema. Keep the bundle outside this code repository. No data are downloaded automatically.

## Requirements

Install Python 3.10+ packages with `python -m pip install -r requirements.txt`. For fitting, install R and the R packages `survival`, `splines`, `cluster`, and `parallel`. `splines` and `parallel` are bundled with R; ensure `survival` and `cluster` are installed. The R executable must be available as `Rscript` or passed with `--rscript`.

## Run

Assuming the requested bundle is unpacked at `/path/to/bundle/input`:

```text
python run.py validate --input /path/to/bundle/input
python run.py fit --input /path/to/bundle/input --output /path/to/new-run --rscript Rscript
```

`--output` must name a new directory outside this repository. The fit command writes `hsc/`, `mpp/`, `figures/`, `report.html`, `input_check.json`, and `run.log` there. It first fits both models and passes **those newly generated result paths** to the renderer. A failed fit stops before rendering. The log contains each invoked command and any R diagnostics.

For a quick rendering check using results that already exist outside the repository:

```text
python run.py render --input /path/to/bundle/input --output /path/to/new-render \
  --hsc-results /path/to/existing-hsc-results --mpp-results /path/to/existing-mpp-results
```

That command labels the report as an external-results rendering check and does **not** claim to have fitted models. The HSC output must contain `model_results/CLP_statistics.csv` and `model_results/LK_statistics.csv`; the MPP output must contain `all_fit_statistics.csv`.

If the supplied data are integrated but still on the raw intensity scale, explicitly prepare the model input:

```text
python code/prepare_input.py --raw-matrix /path/to/raw_intensity.csv.gz \
  --cells /path/to/cells_with_TIC.csv --tic-reference <TIC_REFERENCE> \
  --output /path/to/bundle/input/log_intensity.csv.gz
```

Replace <TIC_REFERENCE> with the documented prefilter TIC median provided with the data. Preserve the matrix's feature IDs and use the same cell IDs throughout the bundle.

## Interpretation

The HSC model treats zero as left-censored on the log1p scale, fits a Gaussian Tobit model with batch and a natural spline of path-normalized pseudotime, and applies BH correction across valid feature × branch tests. The MPP fits use separate within-scope BH families. An invalid fit is recorded separately from a nonsignificant fit. Reported associations are conditional on the supplied trajectory and do not establish branch differences, donor-level replication, causality, or metabolite structural identity.

Keep the exact input bundle, software versions, and run log with any published result.

## Numerical configuration

HSC fitting uses path-specific time normalization and batch weights. Joint clustering standardizes each path separately before combining curves. Numerical and clustering parameters are defined in `code/refit_hsc.R`. Set the optional `HSC_RANDOM_SEED` environment variable to record a seed for stochastic extensions. The current Tobit and PAM calculations are deterministic and do not require a random seed.
