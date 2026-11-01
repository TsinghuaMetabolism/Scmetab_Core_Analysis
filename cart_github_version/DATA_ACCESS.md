# External data and parameter contract

The pre-publication code snapshot contains no study data. Raw and processed
matrices, cell annotations, fitted models, numerical results and the completed
parameter settings will be supplied separately following publication. No data
are downloaded automatically. Keep data and completed settings outside this
code directory.

## Panel checkpoint bundle

`run.py render --input BUNDLE` reads files under
`BUNDLE/figure4/figure4b/OriginalData/`, and analogously for the other panels.
`run.py` lists the exact input filenames in `PANELS`. Tables are tab-separated
unless their extension is `.csv`; JSON files contain model/display metadata.

| Panels | Required data |
|---|---|
| 4b–d | `cells.tsv`: cell_id, round, DDRTree1, DDRTree2, Pseudotime01; `centers.tsv`: center_id, DDRTree1, DDRTree2; `center_edges.tsv`: from, to; `model.json`: root_cell |
| 4c | Also `reference_display.json`: xlim, ylim, root_xy, plot_area_mm |
| 4e | `common_features.tsv`: common_id, label, is_SM, SM_order; `saved_curves.tsv`: common_id, donor, t, fitted_intensity, q |
| 4f | `abundance.tsv`: donor, target_id, native_cell_id, Pseudotime01, abundance; `target_reference.tsv`: reference_row, label; saved spline inputs listed below |
| 4g, 5d | `abundance.tsv`: donor, target_id, cell_id, comparison_group, raw, TIC, label |
| 5e | `cells.tsv`: cell_id, group, actual_day, projection_x_Myr, projection_y_Myr; `display_cells.tsv`: cell_id, plot_x, plot_y; `reference_cells.tsv`: cell_id, DDRTree1, DDRTree2; centers and edges as above; `reference_display.json`: xlim, ylim, root_cell, root_xy, plot_area_mm |
| 5f | `cells.tsv`: unique cell_id, group (Ctr/Myr), t_graph_Myr on the shared frozen [0,1] axis |

For Figure 4f, each `basis/DONOR/` contains `cell_order.tsv` (cell_id,
Pseudotime01), `basis_input.csv` (t), `basis.tsv` (numeric design without
intercept), `basis_grid.csv` (t), `basis_grid.tsv` (prediction design without
intercept), and `basis_metadata.json` (basis, k, penalized). The renderer checks
the saved cubic-regression spline against `mgcv` before fitting. This saved
panel format excludes the intercept; the standalone `basis` command exports
designs including an intercept for the Python fitting functions.

The Figure 4f response is the original-TIC log1p intensity. Figure 4g uses Day0,
Day3 and Day9; Figure 5d uses same-run C0/M0. Target identities and the three
donor lists are supplied in the completed settings, preserving the study's
prespecified five-target display order. Feature IDs and cell IDs must agree
across every file used together. Figure 4e q-values come from each donor's full
dynamic-testing family before target/annotation selection.

## Computational inputs

`CART_analysis.py` accepts cell × feature arrays and aligned metadata.
`tic_log1p` takes original recorded TIC values; `bridge_mean` takes aligned
features and matched source/Bridge anchors. `select_ordering` receives the
prespecified stage pairs, selection threshold and feature count.

`run.py trajectory` reads `ordering.tsv.gz` (features × cells, first column
feature ID) and `cells.tsv` (first column cell ID, numeric round). The matrix
must already be rank-normal transformed, with identical cell order.
`run.py basis` reads `pseudotime.tsv` and `prediction_grid.tsv`, each with a t
column. Dynamic fitting uses the supplied response matrix, feature IDs and
intercept-inclusive spline designs. Fixed-reference projection requires the
reference ranks/detection information, learned loadings/nodes/edges, root,
native reference tree projections and the shared observable feature panel.

The parameter template deliberately uses null fields. Fill parameters from
the approved study configuration; null values are not numerical defaults.
