# Requested data bundle

Expression matrices and cell-level data are not distributed in this code repository. Request the study bundle from the authors. Place its `input/` directory outside the repository and pass that path to `run.py --input`.

```text
input/
  log_intensity.csv.gz
  cells.csv
  annotations.csv
  trajectory/
    final_main_cells.csv
    final_main_edges.csv
    reference_global_combat_umap_cells.csv
  mpp_heatmap/
    cells.csv
    annotated_log_intensity.csv
```

`log_intensity.csv.gz` is a wide CSV: first column `cell_id`, then one column per `feature_id`. Entries are finite, nonnegative **TIC-normalized log1p intensities before ComBat**. Zero means nondetection for these Tobit fits. Do not substitute a ComBat-adjusted matrix, trajectory representation, unnormalized intensity matrix, or instrument RAW data. `cells.csv` must contain unique `cell_id`, `Pseudotime`, `on_CLP_path`, `on_LK_path`, `cell_type`, and `batch`. Every matrix cell must occur exactly once in this file. Path flags identify complete CLP and LK paths, including shared trunk cells.

`annotations.csv` must contain unique `feature_id` and `metabolite` fields and include every main-matrix feature. Annotation names are candidate mass matches, not confirmed structures. The trajectory cell table needs `cell_id`, `cell_type`, `State`, `Pseudotime`, `DDRTree1`, and `DDRTree2`. The edge table needs `x1`, `y1`, `x2`, and `y2`. The UMAP table needs `cell_id`, `cell_type`, `UMAP1`, and `UMAP2`. These coordinate tables must cover exactly the cells in `cells.csv`.

The MPP `cells.csv` needs unique `cell_id`, `batch`, `cell_type` (all `MPP`), `State`, and `Pseudotime`. The MPP matrix must start with `cell_id` followed by annotated feature IDs on the same TIC-log1p scale. Its cell rows must match the MPP metadata in order. The MPP model uses S1, S2, and S3 scopes and excludes a feature from a scope if any batch has no detections. This is an analysis rule, not a missing-data imputation step.

Trajectory and UMAP coordinates are fixed inputs. The expression matrices, annotations, and coordinate tables must use matching identifiers.
