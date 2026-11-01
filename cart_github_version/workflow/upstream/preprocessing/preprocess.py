"""Original-TIC normalization, feature matching and reference ordering."""
import heapq
import numpy as np
import pandas as pd
from scipy import stats


def tic_log1p(raw, tic, scale):
    """Normalize with the original total ion count; preserve zeros and NaNs."""
    raw, tic = np.asarray(raw, float), np.asarray(tic, float)
    if raw.ndim != 2 or tic.shape != (len(raw),):
        raise ValueError('Expected cells x features and one TIC per cell')
    if not np.isfinite(tic).all() or np.any(tic <= 0) or scale <= 0:
        raise ValueError('TIC and scale must be positive and finite')
    if np.isinf(raw).any() or np.any(raw < 0):
        raise ValueError('Raw intensities must be nonnegative; structural gaps use NaN')
    return np.log1p(raw / tic[:, None] * scale)


def bridge_mean(normalized, metadata, anchors):
    """Shift each observable feature to its matching Bridge anchor mean.

    metadata has run and cell_type columns. anchors maps each source run to its
    baseline label; Bridge contains the corresponding baseline labels. Cells and
    columns must already be aligned. Structural gaps remain NaN.
    """
    z = np.asarray(normalized, float)
    if len(z) != len(metadata):
        raise ValueError('Metadata must follow matrix row order')
    y = z.copy()
    shifts = []
    for run, label in anchors.items():
        ii = np.flatnonzero(metadata.run.eq(run))
        si = np.flatnonzero(metadata.run.eq(run) & metadata.cell_type.eq(label))
        ri = np.flatnonzero(metadata.run.eq('Bridge') & metadata.cell_type.eq(label))
        if not len(si) or not len(ri):
            raise ValueError(f'Missing source/Bridge anchor: {run}, {label}')
        cols = np.flatnonzero(np.isfinite(z[ii]).all(0) & np.isfinite(z[ri]).all(0))
        delta = z[np.ix_(ri, cols)].mean(0) - z[np.ix_(si, cols)].mean(0)
        y[np.ix_(ii, cols)] += delta
        shifts.append(pd.DataFrame({'run': run, 'column': cols, 'shift': delta}))
    return y, pd.concat(shifts, ignore_index=True)


def select_ordering(values, stages, feature_ids, pairs, count, alpha):
    """MWU screen: Bonferroni-min P across stage pairs, then feature-wise BH.

    Retain all q < alpha features, then fill from p < alpha in P/ID order.
    The supplied count must accommodate all strict features.
    """
    y, stages = np.asarray(values, float), np.asarray(stages)
    if not np.isfinite(y).all():
        raise ValueError('Ordering uses features observed in every reference cell')
    ps = []
    for first, second in pairs:
        if not np.any(stages == first) or not np.any(stages == second):
            raise ValueError('Each specified stage must contain reference cells')
        ps.append(stats.mannwhitneyu(y[stages == second], y[stages == first],
                  axis=0, alternative='two-sided', method='asymptotic',
                  use_continuity=True).pvalue)
    p = np.minimum(1, len(pairs) * np.min(ps, axis=0))
    table = pd.DataFrame({'feature_id': feature_ids, 'p_mwu': p,
                         'q_mwu': stats.false_discovery_control(p)}).sort_values(['p_mwu', 'feature_id'])
    strict = table.loc[table.q_mwu < alpha, 'feature_id'].tolist()
    candidates = [f for f in table.loc[table.p_mwu < alpha, 'feature_id'] if f not in strict]
    if count < len(strict) or count > len(strict) + len(candidates):
        raise ValueError('Requested ordering count is outside the eligible pool')
    selected = strict + candidates[:count - len(strict)]
    positions = pd.Index(feature_ids).get_indexer(selected)
    transformed = stats.norm.ppf((stats.rankdata(y[:, positions], axis=0) - .5) / len(y))
    return table, selected, transformed


def bounded_run_groups(mz, run_codes, *, ppm):
    """Agglomerate by union span, with a diameter bound and one member per run.

    This is a constrained one-dimensional agglomeration, not unbounded linkage:
    every tentative union is checked against all members through its extrema.
    """
    mz = np.asarray(mz, float)
    codes = np.asarray(run_codes, int)
    if len(mz) == 0:
        return []
    order = np.lexsort((np.arange(len(mz)), codes, mz))
    masses, codes = mz[order], codes[order]
    n = len(masses)
    parent = np.arange(n)
    lo, hi = masses.copy(), masses.copy()
    masks = [1 << int(code) for code in codes]
    queue = []
    for i in range(n):
        j = i + 1
        while j < n:
            distance = (masses[j] - masses[i]) / masses[i] * 1e6
            if distance > ppm:
                break
            if codes[i] != codes[j]:
                queue.append((float(distance), i, j))
            j += 1
    heapq.heapify(queue)

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = int(parent[i])
        return i

    while queue:
        priority, first, second = heapq.heappop(queue)
        a, b = root(first), root(second)
        if a == b or masks[a] & masks[b]:
            continue
        a, b = min(a, b), max(a, b)
        lower, upper = min(lo[a], lo[b]), max(hi[a], hi[b])
        span = (upper - lower) / lower * 1e6
        if span > ppm:
            continue
        if span > priority:
            heapq.heappush(queue, (float(span), a, b))
            continue
        parent[b] = a
        lo[a], hi[a] = lower, upper
        masks[a] |= masks[b]

    groups = {}
    for i in range(n):
        groups.setdefault(root(i), []).append(int(order[i]))
    result = [np.asarray(ids, int) for ids in groups.values()]
    result.sort(key=lambda ids: (float(mz[ids].mean()), int(ids.min())))
    assert sorted(np.concatenate(result).tolist()) == list(range(n))
    for ids in result:
        assert len(set(run_codes[i] for i in ids)) == len(ids)
        assert (mz[ids].max() - mz[ids].min()) / mz[ids].min() * 1e6 <= ppm
    return result
