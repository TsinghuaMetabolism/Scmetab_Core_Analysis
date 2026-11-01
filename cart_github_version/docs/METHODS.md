# CAR-T analysis methods

## Preprocessing and reference geometry

The preprocessing module uses original recorded TIC denominators, a supplied normalization scale and log1p transformation. Original nondetections remain zero; features not measured in a run remain structurally missing. Cross-run matching uses a supplied ppm span bound and at most one feature per run in a matched group. Bridge mean correction shifts each observable feature using the corresponding source/Bridge baseline anchors.

Ordering features are screened with stage-pair Mann–Whitney tests, minimum-P Bonferroni adjustment across pairs and feature-wise BH. The prespecified feature count and threshold are external settings. Rank-normal transformed ordering features define the DDRTree geometry; normalized intensities remain the responses for abundance and dynamic analyses. The R module fits the reference, evaluates biological root candidates and exports coordinates, loadings, nodes and native tree projections.

## Dynamic responses and Figure 4e–f

`workflow/upstream/dynamics/fit_dynamics.py` fits an unpenalized cubic-regression spline with an intercept. HC3 covariance gives pointwise 95% intervals for the mean curve. Flatness tests use a null-imposed Mammen wild bootstrap, followed by BH over all finite, variable features in the reference analysis. The spline basis is constructed with `mgcv` using the supplied basis dimension.

Figure 4e uses frozen per-donor curves and their full-family q-values. Rows must have concordant fitted endpoint directions and the prespecified minimum number of donors below the dynamic-q threshold. Mean pairwise correlation orders rows; prespecified SM rows receive display priority after meeting the same inclusion rule. Row Z scaling is computed separately for each donor/feature curve.

Figure 4f fits the five target responses on the original-TIC scale with saved native spline bases and HC3 intervals. Its renderer verifies the saved basis against `mgcv` before fitting.

## Abundance comparisons and Figure 4g/5d

Abundance statistics use two-sided asymptotic Mann–Whitney U tests with tie and continuity correction. All specified donors, targets and group pairs in a panel form one BH family. Figure 4g compares actual Days 0, 3 and 9; Figure 5d compares same-run Bridge C0/M0 groups. Group means and SEM include zero observations.

## Fixed-reference projection and Figure 5e–f

The projection module reuses the Ctr reference feature transform, learned loadings and nodes. It solves conditional coordinates on the shared observable feature panel and projects them onto the frozen graph, including reference-supported leaf extensions. Ctr conditional coordinates use the same panel as Myr for distribution comparisons.

Figure 5e shows native Ctr coordinates as the background and saved Myr projected coordinates. Figure 5f retains all Ctr/Myr cells on their shared frozen [0,1] graph axis, with reflected Gaussian KDE and a two-sided Mann–Whitney test. Groups are not rescaled separately. The analyses describe sampled cell distributions and associations conditional on the reference trajectory.

## Execution and input boundary

`../scripts/run_workflow.py` orchestrates the panel modules listed in `../workflow/provenance/panels.json`. Individual `plot.py` entry points dispatch the same workflow for their panel. The uploaded scientific functions, plotting code and null configuration template are sufficient to inspect the method; executing the study requires the separate data and completed settings described in [DATA_ACCESS.md](DATA_ACCESS.md).
