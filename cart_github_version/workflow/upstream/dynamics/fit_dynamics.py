"""Continuous-response spline fitting, HC3 intervals and Mammen tests."""
import hashlib
import numpy as np
import pandas as pd
from scipy import stats


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
