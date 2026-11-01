from __future__ import annotations
from pathlib import Path
import sys,contextlib,io
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"archive/P59_integration_20260917/vendor"))
from neuroCombat import neuroCombat
def combat_affine(
    x: np.ndarray, batch: np.ndarray, fit_mask: np.ndarray
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], dict]:
    levels = np.unique(batch[fit_mask])
    if len(levels) < 2:
        raise ValueError("At least two batches are required for an anchor correction")
    valid = np.ones(x.shape[1], dtype=bool)
    for level in levels:
        valid &= np.nanvar(x[fit_mask & (batch == level)], axis=0) > 1e-10

    transforms = {
        level: (np.ones(x.shape[1], dtype=float), np.zeros(x.shape[1], dtype=float))
        for level in levels
    }
    covars = pd.DataFrame({"batch": batch[fit_mask]})
    with contextlib.redirect_stdout(io.StringIO()):
        result = neuroCombat(
            x[fit_mask][:, valid].T,
            covars,
            batch_col="batch",
            eb=True,
            parametric=True,
            mean_only=False,
        )
    estimates = result["estimates"]
    pooled_var = estimates["var.pooled"].ravel()
    center = estimates["stand.mean"][:, 0] + estimates["mod.mean"].mean(axis=1)
    learned_levels = np.asarray(estimates["batches"]).astype(str)
    for level in levels:
        level_index = int(np.flatnonzero(learned_levels == level)[0])
        slope_valid = 1.0 / np.sqrt(estimates["delta.star"][level_index])
        offset_valid = (
            center * (1.0 - slope_valid)
            - estimates["gamma.star"][level_index] * np.sqrt(pooled_var) * slope_valid
        )
        slope, offset = transforms[level]
        slope[valid] = slope_valid
        offset[valid] = offset_valid

    invalid = ~valid
    if invalid.any():
        target = np.mean(
            [x[fit_mask & (batch == level)][:, invalid].mean(axis=0) for level in levels], axis=0
        )
        for level in levels:
            slope, offset = transforms[level]
            offset[invalid] = target - x[fit_mask & (batch == level)][:, invalid].mean(axis=0)

    fitted = np.empty_like(x[fit_mask])
    fit_batch = batch[fit_mask]
    for level in levels:
        slope, offset = transforms[level]
        fitted[fit_batch == level] = x[fit_mask & (batch == level)] * slope + offset
    if not np.allclose(fitted[:, valid].T, result["data"], atol=1e-6, rtol=1e-6):
        raise AssertionError("Transferred affine parameters do not reproduce the ComBat training output")

    details = {
        "batches": levels.tolist(),
        "valid_features": int(valid.sum()),
        "constant_feature_fallback": int(invalid.sum()),
        "positive_slopes": bool(all((transforms[level][0] > 0).all() for level in levels)),
    }
    return transforms, details

def day_number(day: np.ndarray) -> np.ndarray:
    return np.array([float(str(value).replace("D", "")) for value in day], dtype=float)

def fit_dual_anchor(dataset: dict, train: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    x = dataset["X"]
    meta = dataset["meta"]
    batch = meta["batch_short"].to_numpy(str)
    day = meta["day"].to_numpy(str)
    mono = np.isin(batch, ["B1", "B2"])
    if train is None:
        train = np.ones(len(meta), dtype=bool)

    d0_mask = train & (day == "D0")
    d0_transforms, d0_details = combat_affine(x, batch, d0_mask)
    baseline = x.copy()
    for level, (slope, offset) in d0_transforms.items():
        baseline[batch == level] = x[batch == level] * slope + offset

    second_mask = train & mono & (day == dataset["anchor_day"])
    second_transforms, second_details = combat_affine(baseline, batch, second_mask)

    weights = np.zeros(len(meta), dtype=float)
    dnum = day_number(day)
    weights[mono] = np.clip(dnum[mono] / dataset["anchor_day_number"], 0.0, 1.0)
    corrected = baseline.copy()
    for level, (slope, offset) in second_transforms.items():
        ix = batch == level
        fully_corrected = baseline[ix] * slope + offset
        corrected[ix] += weights[ix, None] * (fully_corrected - baseline[ix])

    assert np.isfinite(corrected).all()
    model = {
        "method": "Two-stage state-stratified location-scale neuroCombat",
        "preprocessing": "log1p(raw intensity / recorded TIC * 1e7)",
        "stage_1": "Empirical-Bayes location-scale ComBat fitted on D0 from B1/B2/B3 and transferred to all cells.",
        "stage_2": (
            f"Empirical-Bayes location-scale ComBat fitted on monotherapy {dataset['anchor_day']} from B1/B2. "
            f"Applied only to monotherapy trajectories with weight clip(day/{int(dataset['anchor_day_number'])}, 0, 1)."
        ),
        "combination_batch": "B3 treated cells receive only the D0-derived correction.",
        "D0_fit": d0_details,
        "second_anchor_fit": second_details,
        "anchor_day": dataset["anchor_day"],
    }
    return corrected, model
