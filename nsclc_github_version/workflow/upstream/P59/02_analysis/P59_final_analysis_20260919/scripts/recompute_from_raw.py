"""Recompute TIC+log1p and D0+D3 correction from raw files using frozen QC/feature decisions."""
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
from integration_core import fit_dual_anchor
B=Path(__file__).resolve().parents[1];out=B/'recomputed';out.mkdir(exist_ok=True)
meta=pd.read_csv(B/'data/P59_cells_master.csv');f=pd.read_csv(B/'data/P59_features.csv')
mp=pd.read_csv(B/'tables/integration/feature_mapping.csv')
files={'B1_Osi_D8':'P59-D0-D3-D8-20260909-Annotation.h5ad','B2_Osi_D14':'P59-D0-D3-D14-Annotation-20260827.h5ad','B3_Osi_TCD_D8':'P59-perturb-D0-D3-D8-Annotation-20260823.h5ad'}
x=np.empty((len(meta),len(f)),dtype=float)
for b,name in files.items():
    raw=ad.read_h5ad(B.parents[1]/'01_raw/h5ad'/name);ii=np.where(meta.batch.eq(b))[0]
    rr=raw.obs_names.get_indexer(meta.iloc[ii].original_cell_id.astype(str));assert (rr>=0).all()
    mm=mp[mp.batch.eq(b)].set_index('feature_id').loc[f.matrix_feature_id]
    cc=mm.source_column.astype(int).values
    vals=raw.X[rr][:,cc];vals=vals.toarray() if hasattr(vals,'toarray') else np.asarray(vals)
    # QC did not change feature intensities; recorded TIC, not feature-sum TIC.
    x[ii]=np.log1p(vals/meta.iloc[ii].TIC.to_numpy()[:,None]*1e7)
np.testing.assert_allclose(x,np.load(B/'data/P59_TIC_log1p_before.npy'),atol=1e-10,rtol=1e-10)
y,model=fit_dual_anchor({'X':x,'meta':meta,'anchor_day':'D3','anchor_day_number':3})
np.testing.assert_allclose(y,np.load(B/'data/P59_corrected.npy'),atol=1e-8,rtol=1e-8)
np.save(out/'TIC_log1p_before.npy',x);np.save(out/'corrected.npy',y)
(out/'recompute_validation.json').write_text(json.dumps({'raw_to_normalized_match':True,'refitted_correction_match':True,'normalization_max_abs_error':float(np.max(abs(x-np.load(B/'data/P59_TIC_log1p_before.npy')))),'correction_max_abs_error':float(np.max(abs(y-np.load(B/'data/P59_corrected.npy'))))},indent=2))
print((out/'recompute_validation.json').read_text())
