# Restricted input contract

No biological matrices, cell-level tables, fitted objects, figures or results are included in this code-only directory.

The authors should supply an approved bundle containing:

1. Annotated cell-by-feature intensity matrices, with explicit raw/normalized/logged scale and recorded per-cell TIC.
2. Cell ID, original cell ID, donor, sample, batch, treatment series, actual treatment status, day and independent biological replicate metadata.
3. Stable feature IDs, original source-column mapping, m/z, candidate annotations and their ambiguity.
4. QC decisions and batch-feature matching for exact historical reproduction.
5. Reference ordering/scaling, fixed trajectories, current restricted query projections, current GAM inputs and cross-donor mapping/curve checkpoints where required.
6. Frozen compound/pathway annotation resources, enrichment background and source checksum manifest.

`restricted_files.tsv` lists the data/resource/checkpoint paths in the archived manifest without copying file contents or local source paths. This is a delivery inventory, not a declaration that all listed files are raw data. Source paths must remain relative to the approved bundle root.

P41 B2/B3 original H5AD and original annotation CSV were not present in the supplied archive. No replacement data are fabricated. Authors must provide missing inputs and finish the connected raw-restart validation before claiming that original matrices alone reproduce every panel.

Raw instrument peak extraction and chemical identification are upstream of this repository. Do not interpret a candidate feature annotation as a uniquely identified metabolite.

## External annotation paths

For the archived P41 integration module, set `NSCLC_P41_ANNOTATION_CSV` to the original annotation CSV and `NSCLC_RESOURCE_DIR` to the annotation-resource directory. Defaults are `inputs/original_annotation.csv` and `inputs/annotation_resources` relative to that module root.
