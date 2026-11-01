# Methods and fixed decisions

- Normalize measured intensities using log1p(raw / recorded TIC * 1e7). Preserve structural missingness; do not renormalize corrected values.
- Use biologically shared anchors for empirical-Bayes location/scale correction. P59 uses D0 and reference D3; P39/P41 use endpoint-informed variants. Treatment-specific correction is not unrestricted batch mixing.
- Current P39 reference uses 159 endpoint-consistent ordering features, not the older top100 trajectory. P41 stage-based ordering and P59 frozen top100 inputs are distinct historical profiles.
- Fit reference trajectories before projecting query cells. Preserve the exact current query-coordinate version used for downstream inference.
- Gaussian GAM uses k=6, cubic regression splines and REML in the current workflow. Include/exclude-D0 analyses are separate. Apply BH to the full declared feature family before selecting named metabolites.
- Shared38 selection requires concordant non-Complex shape, all pairwise curve correlations >=0.5 and original GAM q<0.05 in at least two donors. These counts are fixture outcomes, not targets for new data.
- KEGG ORA counts unique compounds. Historical background2901 and the21-hit versus86-eligible pathway BH families are reported separately.
- Distribution plots use fixed-reference scaling and reflected KDE. Pooled groups have different time compositions; same-day comparisons remain necessary.
- A single-State unbranched trajectory can reverse pseudotime on refit. The orientation adapter preserves the native model and tests Monocle reverse=False/True using only D0/latest-day median ordering. Native plot coordinates are not rotated or swapped.

## Reproducibility and interpretation

The supported checkpoint-level operations have been numerically compared with the supplied reference results. This does not validate arbitrary new datasets or a connected raw-to-all-panels restart.

Both inclusion and exclusion of D0 must be reported. Choline association in the combination-treatment group is sensitive to this choice and should not be described as significant across both settings. Obtain the corresponding statistics from each regenerated full-feature result table; do not transfer statistics between model versions.

Within-group trajectory association and between-group curve differences are distinct hypotheses. A significant association in one group and a nonsignificant association in another does not establish a group difference. These are exploratory cell-level results from a single donor for the P59 comparison; cells are not independent biological replicates. Pointwise curve intervals do not propagate uncertainty from upstream correction, trajectory estimation or query mapping. Differences in sampling-day composition limit causal interpretation.

## Before publishing

Authors must specify the data access contact/terms, choose a code license, supply missing P41 raw inputs if a complete restart is claimed, and validate replacement of current fixed projection/shared-curve checkpoints. Do not advertise unsupported execution levels.
