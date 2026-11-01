"""Main-panel abundance tests, cross-donor selection and pseudotime distributions."""
import numpy as np
import pandas as pd
from scipy import stats
from workflow.upstream.preprocessing.preprocess import tic_log1p


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
