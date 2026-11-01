"""Portable redraw of Figure 3 b-i from frozen, bundled numerical inputs.
No upstream fitting or reference-coordinate changes are performed here.
"""
from pathlib import Path
import argparse,json,sys,os,hashlib,subprocess,html
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
P=json.loads((ROOT/'provenance/panels.json').read_text(encoding='utf8'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
COLORS={'D0':'#2166AC','D3':'#2CB1BC','D8':'#F6AD31','D14':'#D1493F'}
def read(letter,name):return pd.read_csv(ROOT/P[letter]/'inputs'/name)
def save(fig,letter,name):
    for ext in ['png','svg','pdf']:fig.savefig(ROOT/P[letter]/'figures'/f'{name}.{ext}',dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def trajectory(letter):
    donor={'b':'P59','c':'P41','d':'P39'}[letter];d=read(letter,'Osi_cells_and_trajectory.csv');e=read(letter,'fixed_reference_tree_edges.csv')
    assert d.cell_id.is_unique and len(d)=={'b':1234,'c':886,'d':1068}[letter]
    fig,ax=plt.subplots(figsize=(5.5,5.7),layout='constrained')
    for r in e.itertuples():ax.plot([r.x0,r.x1],[r.y0,r.y1],color='#647887',lw=.55,zorder=1)
    ax.scatter(d.trajectory_1,d.trajectory_2,c=d.day.map(COLORS),s=15,alpha=.85,edgecolors='#333333',linewidths=.2,zorder=2)
    xy=d[['trajectory_1','trajectory_2']].to_numpy();mid=(xy.max(0)+xy.min(0))/2;span=np.ptp(xy,axis=0).max()*.55
    ax.set(xlim=(mid[0]-span,mid[0]+span),ylim=(mid[1]-span,mid[1]+span),xlabel='Component 1',ylabel='Component 2',title=f'{donor} | Osi reference | n={len(d):,}');ax.set_aspect('equal')
    for day,color in COLORS.items():ax.scatter([],[],c=color,s=20,label=f'{day} (n={sum(d.day==day)})')
    ax.legend(ncol=2,loc='upper center',bbox_to_anchor=(.5,-.14),frameon=False,fontsize=9)
    save(fig,letter,'trajectory');d.to_csv(ROOT/P[letter]/'tables/plot_cells.csv',index=False);e.to_csv(ROOT/P[letter]/'tables/plot_edges.csv',index=False)
def heatmap():
    all_=read('e','shared_features_ordered.csv')
    keep=(all_.P39_shape==all_.P41_shape)&(all_.P39_shape==all_.P59_shape)&(all_.P39_shape!='Complex')
    keep&=all_[['cor_P39_P41','cor_P39_P59','cor_P41_P59']].min(axis=1)>=.5
    keep&=(all_[['P39_q','P41_q','P59_q']]<.05).sum(axis=1)>=2
    d=all_[keep].copy();expected=read('e','final_38_features.csv')
    assert len(d)==38 and set(d.shared_id)==set(expected.shared_id)
    d=d.set_index('shared_id').loc[expected.shared_id].reset_index()
    d.to_csv(ROOT/P['e']/'tables/selected_38_recomputed.csv',index=False)
    for mode in ['all','Up','Down']:
        rows=d if mode=='all' else d[d.P39_shape==mode]
        fig,axs=plt.subplots(1,3,figsize=(12.3,max(6,len(rows)*.23)),sharey=True,layout='constrained')
        for ax,donor in zip(axs,['P39','P41','P59']):
            z=read('e',donor+'_shared_rowZ.csv').set_index('shared_id')
            mat=z.loc[rows.shared_id,[f'T{i:03d}' for i in range(1,161)]].to_numpy(float)
            im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2,interpolation='nearest');ax.set_title(donor);ax.set_xticks([0,159],['Start','End']);ax.set_xlabel('Pseudotime')
            if mode=='all':ax.axhline(16.5,color='white',lw=3)
        labels=[str(x).replace('\t',' ') for x in rows.annotation]
        axs[0].set_yticks(range(len(rows)),labels,fontsize=8)
        fig.colorbar(im,ax=axs,shrink=.3,pad=.02,label='Saved row Z');save(fig,'e','heatmap_'+mode)
def enrichment():
    # The independent source validator recalculates correlations, ORA and both BH families.
    panel=ROOT/P['f']
    subprocess.run([sys.executable,'-X','utf8',str(panel/'scripts/02_validate_results.py')],check=True)
    d=pd.read_csv(panel/'results/tables/at_least_two_donors_plot_data.csv').sort_values('p.adjust')
    fig,ax=plt.subplots(figsize=(10.5,5.8),layout='constrained');x=-np.log10(d['p.adjust'].clip(lower=1e-300));y=np.arange(len(d))
    ax.hlines(y,0,x,color='#cbd6e1');sc=ax.scatter(x,y,s=d.Count*40,c=d.FoldEnrich,cmap='YlOrRd',edgecolors='#8393a0',linewidths=.4)
    ax.set_yticks(y,d.Description);ax.invert_yaxis();ax.axvline(-np.log10(.05),ls='--',lw=.8,color='#af3654');ax.set(xlabel='-log10(BH FDR; hit pathways)',title='Shared 38-feature set | KEGG enrichment')
    fig.colorbar(sc,ax=ax,label='Fold enrichment',shrink=.65);save(fig,'f','enrichment');d.to_csv(panel/'tables/plot_data.csv',index=False)
def curves():
    c=read('g','fitted_curves_95CI_both_scenarios.csv');names=read('g','name_mapping.csv');pal=['#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#56B4E9','#6A51A3']
    for scenario in ['exclude_D0','include_D0']:
        fig,axs=plt.subplots(2,2,figsize=(10.8,8.5),layout='constrained')
        for row,sel in enumerate([names.iloc[2:],names.iloc[:2]]):
            d=c[(c.scenario==scenario)&c.feature_id.isin(sel.feature_id)];lo=d.lower.min();hi=d.upper.max();pad=(hi-lo)*.06
            for ax,condition in zip(axs[row],['Osi','Osi+TCD']):
                for i,r in enumerate(sel.itertuples()):
                    a=d[(d.condition==condition)&(d.feature_id==r.feature_id)].sort_values('t')
                    ax.plot(a.t,a.fit,color=pal[i],label=r.short_name,lw=1.8);ax.fill_between(a.t,a.lower,a.upper,color=pal[i],alpha=.15)
                ax.set(xlim=(0,1),ylim=(lo-pad,hi+pad),xlabel='Pseudotime',ylabel='Relative intensity',title=condition)
                ax.legend(loc='upper center',bbox_to_anchor=(.5,-.2),frameon=False,fontsize=8,ncol=3 if row==0 else 2)
        save(fig,'g','GAM_'+scenario)
def projection():
    d=read('h','trajectory_cells.csv').sort_values('draw_order');e=read('h','fixed_reference_tree_edges.csv')
    fig,ax=plt.subplots(figsize=(6.8,6.7),layout='constrained')
    for r in e.itertuples():ax.plot([r.x0,r.x1],[r.y0,r.y1],c='#647887',lw=.55,zorder=1)
    col=np.where(d.day!='D8','#D9D9D9',np.where(d.cohort=='Osi_reference','#0072B2','#D55E00'))
    ax.scatter(d.trajectory_1,d.trajectory_2,c=col,s=15,alpha=.85,edgecolors='#555555',linewidths=.18,zorder=2);ax.set_aspect('equal')
    ax.set(xlabel='Component 1',ylabel='Component 2',title='P59 | Day8 on the fixed Osi reference')
    for color,label in [('#0072B2','Osi D8 (n=196)'),('#D55E00','Osi+TCD D8 (n=154)'),('#D9D9D9','Other cells (n=1,295)')]:ax.scatter([],[],c=color,s=20,label=label)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False,ncol=2,fontsize=9);save(fig,'h','Day8_projection')
    d.to_csv(ROOT/P['h']/'tables/plot_cells.csv',index=False)
def distribution():
    sys.path.insert(0,str(ROOT/P['i']/'scripts/toolkit'))
    from pseudotime_plot import plot_pseudotime
    d=read('i','trajectory_cells.csv');ref=d[d.cohort=='Osi_reference'].pseudotime
    d['reference_time']=(d.pseudotime-ref.min())/(ref.max()-ref.min());d['condition']=d.cohort.map({'Osi_reference':'Osi','TCDOsi_query':'Osi+TCD'})
    original=Figure.savefig
    def padded(self,*a,**kw):kw.update(bbox_inches='tight',pad_inches=.15);return original(self,*a,**kw)
    Figure.savefig=padded
    validation=[]
    try:
        for scenario in ['all','exclude_D0','D8']:
            sub=d if scenario=='all' else d[d.day!='D0'] if scenario=='exclude_D0' else d[d.day=='D8']
            prefix=ROOT/P['i']/'figures'/('distribution_'+scenario)
            r=plot_pseudotime(sub,prefix,group_col='condition',time_col='reference_time',groups=['Osi','Osi+TCD'],normalize=False,title='P59 | Reference-position distributions | '+scenario)
            for suffix in ['statistics.json','curves.csv','bins.csv','used_cells.csv']:
                src=Path(str(prefix)+'_'+suffix);dest=ROOT/P['i']/'tables'/src.name;dest.write_bytes(src.read_bytes());src.unlink()
            bins=pd.read_csv(ROOT/P['i']/'tables'/f'distribution_{scenario}_bins.csv');cur=pd.read_csv(ROOT/P['i']/'tables'/f'distribution_{scenario}_curves.csv')
            assert np.allclose(bins.groupby('group').percent.sum(),100)
            assert all(abs(np.trapezoid(v.density,v.pseudotime)-1)<1e-9 for _,v in cur.groupby('group'))
            assert r['n_input']==sum(r['n']) and r['excluded_invalid_or_unmapped']==0
            sub.to_csv(ROOT/P['i']/'tables'/f'distribution_{scenario}_cell_ids.csv',index=False)
            validation.append({'scenario':scenario,'n':r['n'],'MWU_p':r['MWU_p'],'stars':r['stars'],'histogram_percent_sums':bins.groupby('group').percent.sum().to_dict(),'density_integral_check':True})
    finally:Figure.savefig=original
    (ROOT/'provenance/distribution_validation.json').write_text(json.dumps(validation,indent=2),encoding='utf8')
def build_html():
    def figure(letter,key,title,caption):
        folder=P[letter];src=f'{folder}/figures/{key}.png'
        svg=f'{folder}/figures/{key}.svg'
        return f'<section class="panel panel-{letter}"><h2><b>{letter}</b> {title}</h2><img src="{src}" alt="{html.escape(title)}" loading="lazy" tabindex="0"><p>{caption}</p><div class="links"><a href="{svg}">SVG</a><a href="{src}">PNG</a><a href="{folder}/README.html">Inputs and code</a></div></section>'
    style='''*{box-sizing:border-box}body{margin:0;background:#f6f8fc;color:#263747;font:16px/1.6 Arial,sans-serif}header,main,footer{max-width:1500px;margin:auto;padding:28px}header{background:#fff;border-bottom:3px solid #8061a2}h1{font-size:32px;color:#0068b7;margin:0}h2{font-size:19px;color:#0068b7;margin:0 0 12px}h2 b{font-size:25px;margin-right:10px;color:#172c3c}.intro{max-width:1100px}.figure-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:20px}.panel{background:#fff;padding:20px;min-width:0;border:1px solid #dbe4ed}.panel img{display:block;width:100%;height:auto;cursor:zoom-in}.panel p{font-size:14px}.panel-e{grid-column:span 2}.panel-f{align-self:start}.panel-g{grid-column:span 2}.panel-h{align-self:start}.panel-i{grid-column:span 2}.methods{background:white;padding:24px;margin-top:25px}a{color:#0068b7}a:focus,img:focus{outline:3px solid #7954a2}.links{display:flex;gap:16px;font-size:13px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border-bottom:1px solid #dce5ee;text-align:left;padding:10px;vertical-align:top}.scroll{overflow:auto}select{padding:9px;margin:0 0 14px;max-width:100%}dialog{width:96vw;height:94vh;border:0;padding:20px}dialog img{display:block;max-width:100%;height:auto;margin:auto}dialog button{position:sticky;top:0;padding:10px;background:#fff;border:1px solid #ccc;cursor:pointer}dialog::backdrop{background:#112435bb}@media(max-width:850px){.figure-grid{grid-template-columns:1fr}.panel-e,.panel-g,.panel-i{grid-column:auto}header,main,footer{padding:16px}h1{font-size:25px}}'''
    s='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>NSCLC Figure 3 | Reproducible analysis package</title><style>'+style+'</style><header><h1>NSCLC: continuous metabolic trajectories and treatment comparison</h1><p class="intro">Figure 3, panels b–i. Three independent donor references, shared metabolic dynamics, and P59 combination-treatment projections. Clinical information and panel j are excluded.</p><p>Final display: native Component 1 / Component 2, equal coordinate units, no rotation or axis swap. Click a figure to enlarge it.</p></header><main><div class="figure-grid">'
    for letter,donor in [('b','P59'),('c','P41'),('d','P39')]:s+=figure(letter,'trajectory',donor+' reference trajectory','Frozen Osi cell coordinates and tree edges. Sampling day is shown by color; the reference was not refitted.')
    s+=figure('e','heatmap_all','Shared metabolic dynamics','38 features: 17 Up and 21 Down. Same non-Complex shape in all donors, all pairwise curve correlations ≥0.5, and original GAM q&lt;0.05 in at least two donors. Saved row Z values; donor-specific axes are not aligned biological time.')
    s+=figure('f','enrichment','KEGG enrichment of the same 38 features','15 unique mapped compounds; 2,901-compound background. Seven pathways pass BH across 21 hit pathways, while two pass BH across all 86 eligible pathways. Point area indicates hit count. Enrichment does not establish pathway activation.')
    s+=figure('g','GAM_exclude_D0','P59 choline-related metabolic curves','Seven PC-related features (top), Choline and GPC (bottom). Frozen group-specific GAM curves with pointwise 95% intervals. D0 is excluded in this primary view; the included-D0 alternative is available below.')
    s+=figure('h','Day8_projection','Day8 treatment projection','Osi n=196; Osi+TCD n=154; other cells are grey. Current restricted projection onto the unchanged Osi reference. Trajectory position alone does not identify DTP cells.')
    s+=figure('i','distribution_all','Pseudotime-position distributions','All cells include D0 and unequal sampling-day compositions. Shared reference scaling, shared bandwidth, reflected KDE at 0/1. Stars use two-sided asymptotic MWU raw p-values. These are exploratory cell-level comparisons.')
    s+='<section class="methods"><h2>Alternative comparisons</h2><p><a href="'+P['i']+'/figures/distribution_D8.png">Day8 only</a>: MWU p=0.663 (ns), despite a KS shape difference.</p><p><a href="'+P['i']+'/figures/distribution_exclude_D0.png">Exclude D0</a>: MWU p=0.0144.</p><p><a href="'+P['g']+'/figures/GAM_include_D0.png">GAM including D0</a>: sensitivity analysis, not interchangeable with the primary view.</p><p><a href="'+P['e']+'/figures/heatmap_Up.svg">17 Up features</a> / <a href="'+P['e']+'/figures/heatmap_Down.svg">21 Down features</a></p></section></div>'
    s+='''<section class="methods"><h2>Reproduction levels and raw-data availability</h2><div class="scroll"><table><tr><th>Level</th><th>Included material</th><th>Status</th></tr><tr><td>Figure replay</td><td>Per-panel inputs, portable plotting code, original source code, SVG/PNG/PDF outputs and statistics</td><td>b–i regenerated in this package</td></tr><tr><td>Original supplied data</td><td>P39: 2 H5AD files; P59: 3 H5AD files; P41: B1 H5AD</td><td>Copied byte-for-byte under upstream/&lt;donor&gt;/01_raw; these are the supplied annotated H5AD files, not instrument acquisition files</td></tr><tr><td>Upstream analysis</td><td>Selected integration, ordering-feature, DDRTree, projection and GAM workflows; frozen matrices, objects, annotation/mapping tables and historical environments</td><td>Archived as source and checkpoints, not newly refitted. P41 B2/B3 original H5AD and original annotation CSV remain missing.</td></tr></table></div><p>The P41 imported snapshot has altered missing-value representation and is not a substitute for an untouched original H5AD. Historical upstream scripts may require path/environment adaptation; use the portable redraw entry below for exact frozen-input figures.</p><p><a href="README.html">Package guide and commands</a> · <a href="provenance/source_manifest.csv">Original-source checksums</a> · <a href="provenance/file_integrity.json">Integrity validation</a> · <a href="provenance/analysis_entrypoints.csv">Upstream entry points</a></p></section></main><footer>NSCLC Figure 3 reproducible package · 2026-10-14 · Local/offline figures and links. Candidate annotations and cell-level associations do not establish causal mechanisms.</footer><dialog id="zoom"><button id="close">Close</button><img id="large" alt="Enlarged scientific figure"></dialog><script>const z=document.getElementById('zoom');document.querySelectorAll('.panel img').forEach(im=>{const open=()=>{document.getElementById('large').src=im.src;document.getElementById('large').alt=im.alt;z.showModal()};im.onclick=open;im.onkeydown=e=>{if(e.key==='Enter')open()}});document.getElementById('close').onclick=()=>z.close();z.onclick=e=>{if(e.target===z)z.close()};</script></html>'''
    (ROOT/'index.html').write_text(s,encoding='utf8')
def verify():
    rows=pd.read_csv(ROOT/'provenance/source_manifest.csv');missing=[];bad=[]
    for r in rows.itertuples():
        p=ROOT/r.package_path
        if not p.exists():missing.append(r.package_path);continue
        h=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(2**20),b''):h.update(block)
        if h.hexdigest()!=r.sha256:bad.append(r.package_path)
    result={'files':len(rows),'missing':missing,'changed':bad,'passed':not missing and not bad}
    (ROOT/'provenance/file_integrity.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    assert result['passed'],result
    print(json.dumps(result))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['redraw','verify','html'],default='redraw',nargs='?');ap.add_argument('--panel',choices=list(P)+['all'],default='all');args=ap.parse_args()
    if args.action=='verify':verify()
    elif args.action=='html':build_html()
    else:
        for letter in list(P) if args.panel=='all' else [args.panel]:
            print('Drawing',letter,flush=True)
            if letter in 'bcd':trajectory(letter)
            else:{'e':heatmap,'f':enrichment,'g':curves,'h':projection,'i':distribution}[letter]()
        build_html()
