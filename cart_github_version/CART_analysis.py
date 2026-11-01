"""Computational methods for the CAR-T main figures.

Study inputs and fitted parameter values are supplied separately. Functions use
cell-level data; distribution tests are conditional on the sampled cells and
frozen reference. See README.md for input conventions and figure coverage.
"""
import heapq
import hashlib
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax
from scipy.stats import norm


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


def fit_dynamics(values, feature_ids, basis, prediction_basis, grid, repeats, seed):
    """Unpenalized cr spline, HC3 intervals and Mammen flatness test.

    Basis matrices include the intercept. BH spans every finite, variable
    feature, before selecting annotations or targets for the main figures.
    """
    y = np.asarray(values, float)
    x, xp = np.asarray(basis, float), np.asarray(prediction_basis, float)
    records, curves = [], []
    for j, feature in enumerate(feature_ids):
        response = y[:, j]
        if not np.isfinite(response).all() or np.ptp(response) == 0:
            continue
        beta, covariance, _, _ = fit_ols_hc3(response, x)
        fitted = xp @ beta
        half = 1.96 * np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', xp, covariance, xp), 0))
        feature_seed = (seed + int(hashlib.sha256(str(feature).encode()).hexdigest()[:8], 16)) % (2**32)
        result = mammen_flat_test(response, x, repeats, feature_seed)
        records.append({'feature_id': feature, **result, 'end_minus_start': fitted[-1] - fitted[0]})
        curves.append(pd.DataFrame({'feature_id': feature, 't': grid, 'fitted_intensity': fitted,
                                    'lower95': fitted - half, 'upper95': fitted + half}))
    if not records:
        raise ValueError('No finite, variable feature available')
    tests = pd.DataFrame(records)
    tests['q_mammen'] = stats.false_discovery_control(tests.p_mammen.to_numpy())
    return tests, pd.concat(curves, ignore_index=True)


def abundance_statistics(table, pairs, scale):
    """Figure 4g/5d: all specified donor x target x pair comparisons in one BH family.

    table: donor, target_id, cell_id, comparison_group, raw, TIC. Raw zeros stay
    in every summary and test. Positive rank-biserial means second group higher.
    """
    if table.duplicated(['donor', 'cell_id', 'target_id']).any():
        raise ValueError('Duplicate cell/target rows')
    if not np.isfinite(table[['raw', 'TIC']]).all().all():
        raise ValueError('Use structurally measured target features')
    data = table.copy()
    data['abundance'] = tic_log1p(data.raw.to_numpy()[:, None], data.TIC.to_numpy(), scale).ravel()
    summaries, tests = [], []
    for (donor, target, group), part in data.groupby(['donor', 'target_id', 'comparison_group'], sort=False):
        v = part.abundance.to_numpy()
        summaries.append(dict(donor=donor, target_id=target, comparison_group=group,
                              n_cells=len(v), mean=v.mean(), median=np.median(v),
                              sd=v.std(ddof=1), sem=v.std(ddof=1)/np.sqrt(len(v)),
                              detection_rate=float((part.raw > 0).mean())))
    for (donor, target), part in data.groupby(['donor', 'target_id'], sort=False):
        for first, second in pairs:
            a = part.loc[part.comparison_group.eq(first), 'abundance'].to_numpy()
            b = part.loc[part.comparison_group.eq(second), 'abundance'].to_numpy()
            if len(a) < 2 or len(b) < 2:
                raise ValueError('Each prespecified group needs at least two cells')
            result = stats.mannwhitneyu(a, b, alternative='two-sided', method='asymptotic', use_continuity=True)
            tests.append(dict(donor=donor, target_id=target, group_first=first, group_second=second,
                              n_first=len(a), n_second=len(b), U_first=result.statistic, p_value=result.pvalue,
                              rank_biserial_second_minus_first=1-2*result.statistic/(len(a)*len(b)),
                              mean_difference_second_minus_first=b.mean()-a.mean(),
                              median_difference_second_minus_first=np.median(b)-np.median(a)))
    tests = pd.DataFrame(tests)
    tests['q_value'] = stats.false_discovery_control(tests.p_value.to_numpy())
    return data, pd.DataFrame(summaries), tests


def joint_heatmap(curves, feature_metadata, donors, alpha, minimum_donors):
    """Figure 4e: concordant endpoint direction plus prespecified dynamic q rule.

    curves: common_id, donor, t, fitted_intensity, q (full-family dynamic BH).
    metadata: common_id, is_SM, SM_order. Identities are matched upstream using
    annotations and source masses. Correlation orders rows; it is not a filter.
    """
    audit, output = [], []
    for row in feature_metadata.to_dict('records'):
        groups = [curves.loc[curves.common_id.eq(row['common_id']) & curves.donor.eq(d)].sort_values('t') for d in donors]
        if not all(len(v) > 1 and v.t.tolist() == groups[0].t.tolist() for v in groups):
            raise ValueError('Donor curves require the same relative-pseudotime grid')
        values = np.array([v.fitted_intensity.to_numpy() for v in groups])
        q = np.array([v.q.iloc[0] for v in groups])
        eligible = np.isfinite(values).all() and np.isfinite(q).all() and (values.std(axis=1) > 0).all()
        row.update(selected=False, group='', mean_pairwise_r=np.nan)
        if eligible:
            delta = values[:, -1] - values[:, 0]
            corr = np.corrcoef(values)[np.triu_indices(len(donors), 1)]
            selected = ((delta > 0).all() or (delta < 0).all()) and (q < alpha).sum() >= minimum_donors
            group = 'Concordant net increase' if (delta > 0).all() else 'Concordant net decrease'
            row.update(selected=bool(selected), group=group if selected else '', mean_pairwise_r=float(corr.mean()))
            if selected:
                z = (values-values.mean(axis=1, keepdims=True))/values.std(axis=1, ddof=1, keepdims=True)
                for i, part in enumerate(groups):
                    output.append(part.assign(curve_z=z[i], display_z=np.clip(z[i], -2, 2)))
        audit.append(row)
    audit = pd.DataFrame(audit)
    rows = audit.loc[audit.selected].copy()
    rows['group_order'] = rows.group.map({'Concordant net increase': 0, 'Concordant net decrease': 1})
    rows['SM_first'] = ~rows.is_SM.astype(bool)
    rows['priority_order'] = np.where(rows.is_SM, rows.SM_order, 0)
    rows = rows.sort_values(['group_order', 'SM_first', 'priority_order', 'mean_pairwise_r', 'common_id'],
                            ascending=[True, True, True, False, True], kind='stable')
    rows['row_order'] = np.arange(len(rows))
    long = pd.concat(output, ignore_index=True) if output else pd.DataFrame(columns=[*curves.columns, 'curve_z', 'display_z'])
    return audit, rows, long.merge(rows[['common_id', 'row_order']], on='common_id', validate='many_to_one')


def pseudotime_distribution(reference, treatment, bandwidth_floor):
    """Figure 5f: shared frozen [0,1] axis, MWU and boundary-reflected KDE.

    Both groups use conditional projection with the same observable feature
    panel and the same reference scaling. All cells are retained.
    """
    a, b = np.asarray(reference, float), np.asarray(treatment, float)
    if any(len(v) < 2 or not np.isfinite(v).all() or np.any((v < 0) | (v > 1)) for v in (a, b)):
        raise ValueError('Expected at least two finite positions per group in [0,1]')
    pooled = np.concatenate([a, b])
    bandwidth = max(float(pooled.std(ddof=1) * len(pooled)**(-.2)), bandwidth_floor)
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError('Bandwidth must be positive')
    grid, edges = np.linspace(0, 1, 1001), np.linspace(0, 1, 21)
    curves, bins, summaries = [], [], []
    for name, values in [('Ctr', a), ('Myr', b)]:
        density = np.zeros(len(grid))
        for start in range(0, len(values), 1000):
            block = values[start:start+1000]
            for centers in (block, -block, 2-block):
                density += np.exp(-.5*((grid[:, None]-centers[None, :])/bandwidth)**2).sum(axis=1)
        density /= len(values)*bandwidth*np.sqrt(2*np.pi)
        density /= np.trapezoid(density, grid)
        q1, median, q3 = np.quantile(values, [.25, .5, .75], method='linear')
        iqr = q3-q1
        summaries.append(dict(group=name, n=len(values), q1=q1, median=median, q3=q3,
                              whisker_low=values[values >= q1-1.5*iqr].min(),
                              whisker_high=values[values <= q3+1.5*iqr].max()))
        curves.append(pd.DataFrame(dict(group=name, t=grid, density=density, percent_per_bin_width=density*.05*100)))
        counts = np.histogram(values, edges)[0]
        bins.append(pd.DataFrame(dict(group=name, left=edges[:-1], right=edges[1:], n_cells=counts, percent=counts/len(values)*100)))
    mw = stats.mannwhitneyu(b, a, alternative='two-sided', method='asymptotic', use_continuity=True)
    result = dict(U_Myr=float(mw.statistic), P=float(mw.pvalue), q=float(mw.pvalue), BH_family=1,
                  median_difference_Myr_minus_Ctr=float(np.median(b)-np.median(a)),
                  rank_biserial_Myr_minus_Ctr=float(2*mw.statistic/(len(a)*len(b))-1), bandwidth=bandwidth)
    return result, pd.DataFrame(summaries), pd.concat(curves, ignore_index=True), pd.concat(bins, ignore_index=True)


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

def fit_ols_hc3(y,x):
    """Continuous signed outcome, OLS mean and HC3 covariance; no count family."""
    y=np.asarray(y,float);x=np.asarray(x,float)
    if len(y)<=x.shape[1]+1 or np.linalg.matrix_rank(x)!=x.shape[1]:raise ValueError('Insufficient cells or rank-deficient spline model')
    inv=np.linalg.inv(x.T@x);a=inv@x.T;beta=a@y;resid=y-x@beta
    h=np.einsum('ij,ji->i',x,a)
    if np.max(h)>=1-1e-8:raise ValueError('Unstable model leverage')
    u=resid/(1-h);cov=(a*u)@(a*u).T
    return beta,cov,resid,h

def mammen_flat_test(y,x,repeats,seed,chunk_size=512):
    """Null-imposed wild bootstrap, fixed design and unpenalized partial-F statistic.

    y*=mean(y)+HC2-adjusted full-model residual*Mammen weight. Weights have mean
    0, variance 1. Both null and full RSS are recalculated per bootstrap sample.
    This tests all varying coefficients jointly, not significance at every t.
    """
    if repeats<1:raise ValueError('repeats must be positive')
    beta,_,resid,h=fit_ols_hc3(y,x);q,_=np.linalg.qr(x,mode='reduced');n,p=x.shape
    def statistic(values):
        z=q.T@values;ss=np.sum(values**2,axis=0)
        varying=np.maximum(np.sum(z[1:]**2,axis=0),0.)
        rss=np.maximum(ss-np.sum(z**2,axis=0),np.finfo(float).eps)
        return varying/(p-1)/(rss/(n-p))
    # Centering is allowed because the model includes an intercept.
    observed=float(statistic(np.asarray(y)-np.mean(y)))
    adjusted=resid/np.sqrt(1-h);rng=np.random.default_rng(seed)
    sqrt5=np.sqrt(5);low=(1-sqrt5)/2;high=(1+sqrt5)/2;prob_low=(sqrt5+1)/(2*sqrt5)
    exceed=0
    for start in range(0,repeats,chunk_size):
        m=min(chunk_size,repeats-start)
        weights=np.where(rng.random((n,m))<prob_low,low,high)
        # Null intercept cancels from both RSS; omit it for numerical stability.
        star=adjusted[:,None]*weights
        exceed+=int(np.sum(statistic(star)>=observed))
    return {'statistic':observed,'p_mammen':(exceed+1)/(repeats+1),'bootstrap_exceedances':exceed,'bootstrap_repeats':repeats,'p_resolution':1/(repeats+1)}

def reference_transform(values, detected, positive, nzero, nref):
    """Apply the existing detection-aware Ctr midrank transform to each group."""
    z = np.empty_like(values, dtype=float)
    low = np.zeros_like(detected, dtype=bool)
    high = np.zeros_like(detected, dtype=bool)
    for j, ref in enumerate(positive):
        seen = detected[:, j]
        u = np.full(len(values), max(nzero[j] / (2 * nref), .5 / nref))
        x = values[seen, j]
        left = np.searchsorted(ref, x, side='left')
        right = np.searchsorted(ref, x, side='right')
        ranks = np.clip((left + right) / 2, .5, len(ref) - .5)
        u[seen] = (nzero[j] + ranks) / nref
        low[seen, j] = x < ref[0]
        high[seen, j] = x > ref[-1]
        z[:, j] = norm.ppf(np.clip(u, .5 / nref, 1 - .5 / nref))
    return z, low, high

def reduced_objective(z, v, gram, nodes, gamma, sigma):
    """Profile r out; omit the constant x'x when comparing candidate positions."""
    dist2 = np.sum((z[..., None, :] - nodes) ** 2, axis=-1)
    return (np.einsum('...i,ij,...j->...', z, gram, z)
            - 2 * np.sum(z * v, axis=-1)
            - gamma * sigma * logsumexp(-dist2 / sigma, axis=-1))

def value_grad(z, v, gram, nodes, gamma, sigma):
    dist2 = np.sum((nodes - z) ** 2, axis=1)
    r = softmax(-dist2 / sigma)
    f = z @ gram @ z - 2 * z @ v - gamma * sigma * logsumexp(-dist2 / sigma)
    grad = 2 * (gram @ z - v + gamma * (z - r @ nodes))
    return f, grad

def solve_coordinates(x, w, nodes, edge_indices, gamma, sigma):
    """Multi-start alternating minimization, followed by gradient checks."""
    gram = w.T @ w
    inv_a = np.linalg.inv(gram + gamma * np.eye(2))
    v = x @ w
    mapped_nodes = (v[:, None, :] + gamma * nodes[None, :, :]) @ inv_a
    edge_mid = (nodes[edge_indices[:, 0]] + nodes[edge_indices[:, 1]]) / 2
    linear = v @ np.linalg.inv(gram)
    seeds = np.concatenate([
        mapped_nodes, np.broadcast_to(edge_mid, (len(x), len(edge_mid), 2)),
        linear[:, None, :]], axis=1)
    initial_values = reduced_objective(seeds, v[:, None, :], gram, nodes, gamma, sigma)
    z = seeds.copy()
    # Coordinate descent is monotone in the full conditional objective.
    for iteration in range(1000):
        r = softmax(-np.sum((z[:, :, None, :] - nodes) ** 2, axis=-1) / sigma, axis=-1)
        new_z = (v[:, None, :] + gamma * (r @ nodes)) @ inv_a
        step = np.linalg.norm(new_z - z, axis=-1)
        z = new_z
        if step.max() < 1e-11:
            break
    # Only repair starts that have not settled, using the same exact objective.
    repairs = np.argwhere(step >= 1e-8)
    for i, j in repairs:
        opt = minimize(value_grad, z[i, j], args=(v[i], gram, nodes, gamma, sigma), jac=True,
                       method='BFGS', options={'gtol': 1e-9, 'maxiter': 1000})
        z[i, j] = opt.x
    values = reduced_objective(z, v[:, None, :], gram, nodes, gamma, sigma)
    assert np.max(values - initial_values) < 1e-7
    best = values.argmin(axis=1)
    row = np.arange(len(x))
    winner = z[row, best]
    f = values[row, best]
    far = np.linalg.norm(z - winner[:, None, :], axis=-1) > 1e-5
    distinct_values = np.where(far, values, np.inf)
    second = distinct_values.argmin(axis=1)
    second_value = distinct_values[row, second]
    alternative = z[row, second].copy()
    alternative[~np.isfinite(second_value)] = np.nan
    r = softmax(-np.sum((winner[:, None, :] - nodes) ** 2, axis=-1) / sigma, axis=-1)
    mean = r @ nodes
    grad = 2 * (winner @ gram - v + gamma * (winner - mean))
    centered = nodes[None, :, :] - mean[:, None, :]
    covariance = np.einsum('nk,nki,nkj->nij', r, centered, centered)
    hessian = 2 * (gram + gamma * np.eye(2))[None] - 4 * gamma / sigma * covariance
    hessmin = np.linalg.eigvalsh(hessian)[:, 0]
    assert np.max(np.linalg.norm(grad, axis=1)) < 1e-5
    assert np.min(hessmin) > -1e-5
    # Exact bound from log-sum-exp, reported without claiming global optimality.
    node_hard_min = np.min(
        np.einsum('nki,ij,nkj->nk', mapped_nodes, gram, mapped_nodes)
        - 2 * np.sum(mapped_nodes * v[:, None, :], axis=-1)
        + gamma * np.sum((mapped_nodes - nodes[None]) ** 2, axis=-1), axis=1)
    lower_bound = node_hard_min - gamma * sigma * np.log(len(nodes))
    assert np.all(f >= lower_bound - 1e-7)
    return winner, pd.DataFrame({
        'objective_without_x2': f,
        'gradient_norm': np.linalg.norm(grad, axis=1),
        'hessian_min_eigenvalue': hessmin,
        'alternate_objective_gap': np.where(np.isfinite(second_value), second_value - f, np.nan),
        'alternative_Z1': alternative[:, 0], 'alternative_Z2': alternative[:, 1],
        'objective_lower_bound_gap': f - lower_bound,
        'max_node_weight': r.max(axis=1),
        'multistart_count': z.shape[1],
        'max_alternating_iterations': iteration + 1,
        'starts_repaired_with_bfgs': len(repairs),
    })

def freeze_tree(nodes_frame, edges, root_xy, native_projection):
    coordinates={str(k):v for k,v in zip(nodes_frame.index,nodes_frame.to_numpy(float))}
    adjacency={k:[] for k in coordinates}
    for a,b in edges[['from','to']].itertuples(index=False,name=None):
        adjacency[a].append(b);adjacency[b].append(a)
    assert len(edges)==len(coordinates)-1
    leaves=sorted(k for k,v in adjacency.items() if len(v)==1)
    original_branches=sorted(k for k,v in adjacency.items() if len(v)>2)
    root_leaf=min(leaves,key=lambda k:np.linalg.norm(coordinates[k]-root_xy))
    extensions=[]
    root_node=root_leaf
    for leaf in leaves:
        neighbor=adjacency[leaf][0]
        outward=coordinates[leaf]-coordinates[neighbor]
        outward=outward/np.linalg.norm(outward)
        native=native_projection.loc[native_projection.nearest_center.eq(leaf),['P1','P2']].to_numpy()
        extension=max(0.,float(np.max((native-coordinates[leaf])@outward))) if len(native) else 0.
        bound=''
        if extension>1e-12:
            bound=f'Ctr_bound_{leaf}'
            coordinates[bound]=coordinates[leaf]+extension*outward
            adjacency[bound]=[leaf];adjacency[leaf].append(bound)
            if leaf==root_leaf:root_node=bound
        extensions.append(dict(learned_leaf=leaf,bound_node=bound,extension_length=extension,native_cells=len(native),root_side=leaf==root_leaf))
    # Each segment ends at a fork or a leaf. Numbering is deterministic from the root.
    seen=set();distance={root_node:0.};parent={root_node:None};rows=[];segment_count=0
    def walk(node,segment):
        nonlocal segment_count
        assert node not in seen
        seen.add(node)
        children=sorted(k for k in adjacency[node] if k!=parent[node])
        for child in children:
            next_segment=segment
            if segment is None or len(children)>1:
                segment_count+=1;next_segment=f'S{segment_count}'
            parent[child]=node
            length=float(np.linalg.norm(coordinates[child]-coordinates[node]))
            assert length>0
            distance[child]=distance[node]+length
            rows.append(dict(graph_edge=len(rows),from_node=node,to_node=child,
                x=coordinates[node][0],y=coordinates[node][1],xend=coordinates[child][0],yend=coordinates[child][1],
                length=length,root_distance_from=distance[node],root_distance_to=distance[child],tree_segment=next_segment,
                terminal_edge=len(adjacency[child])==1))
            walk(child,next_segment)
    walk(root_node,None)
    assert seen==set(coordinates)
    edge_frame=pd.DataFrame(rows)
    max_depth=max(distance.values())
    node_frame=pd.DataFrame([dict(node_id=k,x=v[0],y=v[1],root_distance=distance[k],t_graph=distance[k]/max_depth,
        degree=len(adjacency[k]),root_node=k==root_node,learned_node=k in nodes_frame.index) for k,v in coordinates.items()])
    return dict(nodes=node_frame,edges=edge_frame,extensions=pd.DataFrame(extensions),
        root_node=root_node,root_leaf=root_leaf,learned_branch_nodes=original_branches,
        learned_leaf_nodes=leaves,max_root_distance=max_depth,total_tree_length=edge_frame.length.sum())

def project_graph(z,tree):
    edges=tree['edges'];a=edges[['x','y']].to_numpy();b=edges[['xend','yend']].to_numpy()
    delta=b-a;lengths=edges.length.to_numpy()
    fraction=np.clip(np.einsum('nej,ej->ne',z[:,None,:]-a,delta)/lengths[None,:]**2,0,1)
    q=a[None]+fraction[...,None]*delta[None]
    dist2=np.sum((z[:,None]-q)**2,axis=2)
    edge=dist2.argmin(axis=1);row=np.arange(len(z));selected=q[row,edge];f=fraction[row,edge]
    distance=edges.root_distance_from.to_numpy()[edge]+f*lengths[edge]
    terminal=(f>=1-1e-10)&edges.terminal_edge.to_numpy()[edge]
    return pd.DataFrame(dict(graph_x=selected[:,0],graph_y=selected[:,1],graph_edge=edge,
        graph_from=edges.from_node.to_numpy()[edge],graph_to=edges.to_node.to_numpy()[edge],
        tree_segment=edges.tree_segment.to_numpy()[edge],graph_edge_fraction=f,
        root_distance=distance,t_graph=np.clip(distance/tree['max_root_distance'],0,1),
        distance_to_tree=np.sqrt(dist2[row,edge]),at_early_endpoint=distance<=1e-10,
        at_late_endpoint=terminal,terminal_endpoint=np.where(terminal,edges.to_node.to_numpy()[edge],'')))
