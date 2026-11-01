"""Isolated historical reproduction with machine-readable evidence; no shared-state writes."""
from pathlib import Path
import argparse, contextlib, csv, hashlib, json, os, platform, runpy, shutil, subprocess, sys, time, traceback

P59 = Path('upstream/P59/02_analysis/P59_final_analysis_20260919')
GAM = Path('upstream/P59/02_analysis/P59_choline_25_GAM_full1928_20261005')

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()

def write(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')

def verify(root, staged=False):
    manifest = root / 'provenance/source_manifest.csv'
    if not manifest.is_file(): raise ValueError('Missing source_manifest.csv')
    rows = list(csv.DictReader(manifest.open(encoding='utf-8-sig')))
    if not rows: raise ValueError('Empty source manifest')
    failures = []; checked=0
    for row in rows:
        p = (root / row['package_path']).resolve()
        if not p.is_relative_to(root): raise ValueError('Manifest path escapes package')
        rel=p.relative_to(root)
        # These are regenerated outputs; source originals were checked separately.
        if staged and ((rel.parts[0].startswith('panel_') and len(rel.parts)>1 and rel.parts[1] in ('tables','figures'))
                       or p.is_relative_to(root/GAM/'tables')): continue
        checked+=1
        if not p.is_file() or digest(p) != row['sha256']: failures.append(row['package_path'])
    if failures: raise ValueError('Input integrity failed: ' + repr(failures[:20]))
    return {'files_checked': checked, 'passed': True, 'generated_outputs_excluded':staged}

def compare_csv(expected, actual, keys, atol=1e-9, rtol=1e-7):
    import numpy as np
    import pandas as pd
    a, b = pd.read_csv(expected), pd.read_csv(actual)
    if set(a.columns) != set(b.columns): raise ValueError(f'Column mismatch: {actual.name}')
    if 't' in keys:
        # Floating prediction-grid coordinates are measurements, not exact string IDs.
        # Match each group's ordered grid positions, then compare t numerically too.
        group=[k for k in keys if k!='t']
        a=a.sort_values(keys).reset_index(drop=True); b=b.sort_values(keys).reset_index(drop=True)
        a['_grid_index']=a.groupby(group,dropna=False).cumcount()
        b['_grid_index']=b.groupby(group,dropna=False).cumcount()
        keys=group+['_grid_index']
    if keys:
        if a.duplicated(keys).any() or b.duplicated(keys).any(): raise ValueError(f'Duplicate comparison keys: {actual.name}')
        a, b = a.set_index(keys).sort_index(), b.set_index(keys).sort_index()
        if not a.index.equals(b.index): raise ValueError(f'Row identity mismatch: {actual.name}')
    if a.shape != b.shape: raise ValueError(f'Shape mismatch: {actual.name}')
    maxima = {}
    for col in a.columns:
        x, y = a[col], b[col]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            np.testing.assert_allclose(x.to_numpy(float), y.to_numpy(float), atol=atol, rtol=rtol, equal_nan=True, err_msg=f'{actual.name}:{col}')
            delta = np.abs(x.to_numpy(float) - y.to_numpy(float))
            maxima[col] = float(np.nanmax(delta)) if np.isfinite(delta).any() else 0.0
        elif not x.fillna('<NA>').astype(str).equals(y.fillna('<NA>').astype(str)):
            raise ValueError(f'Label mismatch: {actual.name}:{col}')
    return {'file': actual.name, 'rows': len(a), 'atol': atol, 'rtol': rtol, 'max_abs_difference_by_column': maxima, 'passed': True}

def run(argv, cwd, log):
    print('Executing', Path(argv[0]).name, Path(argv[1]).name, flush=True)
    if Path(argv[0]).resolve()==Path(sys.executable).resolve() and Path(argv[1]).suffix=='.py':
        # Windows nested Python processes can resolve conflicting BLAS DLLs.
        # Execute the unchanged Python entry in this already-working interpreter.
        old_argv,old_path,old_cwd=sys.argv[:],sys.path[:],Path.cwd()
        try:
            sys.argv=[str(x) for x in argv[1:]]
            sys.path.insert(0,str(Path(argv[1]).parent));os.chdir(cwd)
            with log.open('w',encoding='utf-8') as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                try: runpy.run_path(str(argv[1]),run_name='__main__')
                except SystemExit as e:
                    if e.code not in (0,None): raise RuntimeError(f'Python entry failed ({e.code}); inspect {log}') from e
        finally:
            sys.argv,sys.path=old_argv,old_path;os.chdir(old_cwd)
        return
    env = os.environ.copy(); env['PYTHONIOENCODING'] = 'utf-8'
    # Resolve DLLs from the selected runtime, not another active Conda environment.
    exe=Path(argv[0]).resolve()
    if os.name=='nt':
        prefix=next((p for p in [exe.parent,*exe.parents] if (p/'conda-meta').is_dir()),None)
        if prefix:
            env['PATH']=os.pathsep.join(str(p) for p in [prefix/'Library/bin',prefix,prefix/'Scripts'])+os.pathsep+env.get('PATH','')
            env['CONDA_PREFIX']=str(prefix)
            if exe.name.lower()=='rscript.exe' and (prefix/'Lib/R').is_dir(): env['R_HOME']=str(prefix/'Lib/R')
    with log.open('w', encoding='utf-8') as stream:
        p = subprocess.run([str(x) for x in argv], cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if p.returncode: raise RuntimeError(f'Command failed ({p.returncode}); inspect {log}')

def stage(src, out):
    if src == out or out.is_relative_to(src) or src.is_relative_to(out):
        raise ValueError('Output and package paths must be disjoint, not nested')
    ident = {'source_path': str(src), 'source_manifest_sha256': digest(src/'provenance/source_manifest.csv'),
             'plot_code_sha256': digest(src/'scripts/rebuild.py')}
    marker = out/'run_identity.json'
    if marker.exists():
        if json.loads(marker.read_text(encoding='utf-8')) != ident: raise ValueError('Run identity differs; choose a new output directory')
        if not (out/'package/run.py').is_file(): raise ValueError('Incomplete staged package; choose a new output directory')
    else:
        if out.exists() and any(out.iterdir()): raise ValueError('Output is not empty and is not a managed run')
        out.mkdir(parents=True, exist_ok=True)
        print('Staging isolated package copy', flush=True)
        shutil.copytree(src, out/'package', ignore=shutil.ignore_patterns('__pycache__', '.matplotlib'))
        write(marker, ident)
    return out/'package'

def audit(pkg):
    counts = {d: len(list((pkg/'upstream'/d/'01_raw').rglob('*.h5ad'))) for d in ('P39','P41','P59')}
    return {'raw_h5ad_counts': counts, 'full_three_donor_raw_restart_ready': False,
            'gaps': ['P41 original B2/B3 H5AD and original annotation CSV absent in supplied bundle',
                     'Current restricted query projection and cross-donor curves are checkpoints in the runner',
                     'P59 raw recomputation uses frozen QC and feature mapping'],
            'scope': 'Input integrity and inventory; not scientific refitting'}

def redraw(src, pkg, reports):
    panels = json.loads((pkg/'provenance/panels.json').read_text(encoding='utf-8'))
    # Remove only expected generated main figures in the isolated copy, so stale artifacts cannot pass.
    names = {'b':'trajectory','c':'trajectory','d':'trajectory','e':'heatmap_all','f':'enrichment',
             'g':'GAM_exclude_D0','h':'Day8_projection','i':'distribution_all'}
    for key, name in names.items():
        for ext in ('png','svg'): (pkg/panels[key]/'figures'/f'{name}.{ext}').unlink(missing_ok=True)
    run([sys.executable, pkg/'scripts/rebuild.py', 'redraw'], pkg, reports/'redraw.log')
    checks = []
    for key, name in names.items():
        for ext in ('png','svg'):
            p = pkg/panels[key]/'figures'/f'{name}.{ext}'
            if not p.is_file() or p.stat().st_size < 100: raise ValueError(f'Missing generated figure {p}')
        # Compare every generated CSV for which a baseline exists, including scientific validators.
        for p in sorted((pkg/panels[key]/'tables').glob('*.csv')):
            baseline = src/p.relative_to(pkg)
            if baseline.exists(): checks.append(compare_csv(baseline, p, []))
    from html.parser import HTMLParser
    class Images(HTMLParser):
        def __init__(self): super().__init__(); self.paths=[]
        def handle_starttag(self, tag, attrs):
            if tag=='img':
                val=dict(attrs).get('src')
                if val: self.paths.append(val)
    page = pkg/'index.html'; parser=Images(); parser.feed(page.read_text(encoding='utf-8'))
    if len(parser.paths) != 8 or any(not (pkg/p).is_file() for p in parser.paths): raise ValueError('HTML image integrity failed')
    return {'scope':'Frozen-input b-i reconstruction; no raw-to-figure refit', 'html':str(page),
            'panels':list(names), 'numeric_comparisons':checks,
            'distribution_validation':json.loads((pkg/'provenance/distribution_validation.json').read_text()),
            'html_image_count':len(parser.paths)}

def p59raw(pkg, reports):
    run([sys.executable, pkg/P59/'scripts/recompute_from_raw.py'], pkg/P59, reports/'p59-raw.log')
    result=json.loads((pkg/P59/'recomputed/recompute_validation.json').read_text())
    if not result['raw_to_normalized_match'] or not result['refitted_correction_match']: raise ValueError('Raw numerical validation failed')
    return {'scope':'Original H5AD to normalized and refitted corrected matrices; frozen QC/feature mapping', **result}

def p59gam(src, pkg, reports, rscript):
    run([rscript, pkg/GAM/'scripts/01_fit_full_GAM.R'], pkg/GAM, reports/'p59-gam.log')
    return validate_gam(src,pkg)

def validate_gam(src,pkg):
    checks=[]
    for scenario in ('include_D0','exclude_D0'):
        f=Path('tables')/f'all_1928_GAM_{scenario}.csv'
        checks.append(compare_csv(src/GAM/f, pkg/GAM/f, ['feature_id','scenario'], atol=1e-7, rtol=1e-5))
    f=Path('tables/fitted_curves_95CI_both_scenarios.csv')
    checks.append(compare_csv(src/GAM/f, pkg/GAM/f, ['feature_id','scenario','condition','t'], atol=1e-7, rtol=1e-5))
    return {'scope':'Full-feature GAM refit from saved corrected matrix and current query-coordinate checkpoint', 'comparisons':checks}

def p59trajectory(src, pkg, reports, rscript):
    import numpy as np
    import pandas as pd
    from scipy.spatial import procrustes
    from scipy.stats import spearmanr
    run([rscript, pkg/P59/'scripts/refit_final_monocle2.R', pkg/P59], pkg/P59, reports/'p59-trajectory.log')
    fn=pkg/P59/'refitted_trajectory/P59_pooled_top100_nc50_sig0.1_s20260918_trajectory_cells.csv'
    new=pd.read_csv(fn).set_index('cell_id')
    panel=json.loads((src/'provenance/panels.json').read_text())['b']
    ref=pd.read_csv(src/panel/'inputs/Osi_cells_and_trajectory.csv').set_index('cell_id')
    if set(new.index)!=set(ref.index): raise ValueError('Trajectory cell IDs differ')
    new=new.loc[ref.index]
    _, _, disparity=procrustes(ref[['trajectory_1','trajectory_2']].to_numpy(), new[['DDRTree1','DDRTree2']].to_numpy())
    rho=float(spearmanr(ref.pseudotime, new.Pseudotime).statistic)
    run([rscript,Path(__file__).with_name('orient_linear_reference.R'),pkg/P59],pkg/P59,reports/'p59-orientation.log')
    folder=pkg/P59/'refitted_trajectory'
    orientation=pd.read_csv(folder/'orientation_audit.csv')
    applied=bool(orientation.applied.all())
    oriented_rho=rho;max_delta=None
    if applied:
        oriented=pd.read_csv(folder/'biologically_oriented_trajectory_cells.csv').set_index('cell_id').loc[ref.index]
        oriented_rho=float(spearmanr(ref.pseudotime,oriented.Pseudotime).statistic)
        max_delta=float(np.max(np.abs(ref.pseudotime.to_numpy()-oriented.Pseudotime.to_numpy())))
    same=bool(oriented_rho>=0.99 and disparity<=0.01)
    return {'scope':'DDRTree refit from saved selected/scaled100-feature matrix; no fresh ordering selection or query refit',
            'n_cells':len(new), 'native_pseudotime_spearman':rho, 'pseudotime_spearman':oriented_rho,
            'orientation_applied':applied,'orientation_audit':orientation.to_dict(orient='records'),
            'oriented_pseudotime_max_abs_difference':max_delta,'procrustes_disparity':float(disparity),
            'agreement_rule':{'minimum_spearman':0.99,'maximum_disparity':0.01},
            'agreement_passed':same, 'status':'passed' if same else 'needs_review',
            'note':'Geometry alignment is diagnostic only. Single-State orientation uses experimental D0/endpoint medians, not the historical reference. Native fit retained; frozen replay plots remain unchanged.'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['audit','redraw','p59-raw','p59-gam','p59-trajectory'])
    ap.add_argument('--package',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--rscript',default='Rscript')
    args=ap.parse_args(); src=args.package.resolve(); out=args.output.resolve()
    # Path safety checked before any report or directory write.
    if src==out or out.is_relative_to(src) or src.is_relative_to(out): raise ValueError('Source/output must be disjoint')
    start=time.time(); report={'action':args.action,'status':'running','source':str(src),'output':str(out),
        'python':sys.version,'platform':platform.platform(),'runner_sha256':digest(Path(__file__).resolve()),
        'skill_script_hashes':{p.name:digest(p) for p in Path(__file__).parent.iterdir() if p.is_file()}}
    try:
        integrity=verify(src); pkg=stage(src,out); reports=out/'reports';reports.mkdir(exist_ok=True)
        report['source_integrity']=integrity
        report['staged_integrity']=verify(pkg, staged=True)
        if args.action=='audit': evidence=audit(pkg)
        elif args.action=='redraw': evidence=redraw(src,pkg,reports)
        elif args.action=='p59-raw': evidence=p59raw(pkg,reports)
        elif args.action=='p59-gam': evidence=p59gam(src,pkg,reports,args.rscript)
        else: evidence=p59trajectory(src,pkg,reports,args.rscript)
        report.update(evidence);report['status']=evidence.get('status','passed')
    except Exception as e:
        report.update(status='failed',error=str(e),traceback=traceback.format_exc())
        # Do not write into an unrelated nonempty directory rejected by stage().
        if not (out/'run_identity.json').exists(): print(json.dumps(report,ensure_ascii=False,indent=2));return 1
    report['elapsed_seconds']=round(time.time()-start,2)
    write(out/'reports'/f'{args.action}.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('numeric_comparisons','comparisons','traceback')},ensure_ascii=False,indent=2))
    return 0 if report['status']=='passed' else 2

if __name__=='__main__': raise SystemExit(main())
