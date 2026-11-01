"""Two-group boxplots and reflected-KDE pseudotime distributions."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def significance_stars(p):
    """Two-sided MWU raw P: *** < .001; ** < .01; * < .05; otherwise ns."""
    return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns'


def plot_pseudotime(data, output_prefix, group_col=None, time_col=None,
                    groups=None, group_map=None, normalize=False,
                    title='Pseudotime distributions', colors=('#327FA8','#B96580'),
                    bin_width=.05, bandwidth=None):
    """Read a DataFrame, CSV/TSV, or AnnData/.h5ad obs; export JPG/PNG/SVG and tables.

    Exactly two groups are compared. `groups=[reference, comparison]` selects/order
    groups explicitly. `group_map` merges identities, e.g. C0/C1/C3 -> Ctr.
    Missing group/nonfinite time rows are excluded and counted. Values must be in
    [0,1]; normalize=True applies ONE pooled min-max transform, never per group.
    Returned dict contains paths and statistics. Does not fit a trajectory.
    """
    if isinstance(data,(str,Path)):
        path=Path(data)
        if path.suffix.lower()=='.h5ad':
            import anndata
            obj=anndata.read_h5ad(path,backed='r');df=obj.obs.copy();obj.file.close()
        else:df=pd.read_csv(path,sep='\t' if path.suffix.lower() in ('.tsv','.txt') else ',')
    elif isinstance(data,pd.DataFrame):df=data.copy()
    elif hasattr(data,'obs'):df=data.obs.copy()
    else:raise TypeError('data must be a DataFrame, AnnData, or CSV/TSV/h5ad path')
    def choose(explicit,candidates):
        if explicit:
            if explicit not in df:raise ValueError(f'Column not found: {explicit}')
            return explicit
        for col in candidates:
            if col in df:return col
        raise ValueError(f'Please specify a column; tried {candidates}')
    group_col=choose(group_col,['cell_type','batch','group'])
    time_col=choose(time_col,['pseudotime','Pseudotime','t_graph','Pseudotime01','pseudotime_01'])
    n_input=len(df);g=df[group_col]
    if group_map is not None:g=g.map(group_map)
    t=pd.to_numeric(df[time_col],errors='coerce')
    valid=g.notna() & np.isfinite(t)
    x=pd.DataFrame({'group':g[valid].astype(str),'pseudotime':t[valid].astype(float)})
    excluded_invalid=int((~valid).sum())
    groups=list(map(str,groups)) if groups is not None else list(x.group.unique())
    if len(groups)!=2 or len(set(groups))!=2:raise ValueError('Exactly two groups required; pass groups=[reference,comparison] or group_map to merge labels.')
    excluded_other=int((~x.group.isin(groups)).sum());x=x[x.group.isin(groups)].copy()
    if any((x.group==g).sum()<2 for g in groups):raise ValueError('Each group requires at least two finite observations.')
    lo,hi=float(x.pseudotime.min()),float(x.pseudotime.max())
    if normalize:
        if hi<=lo:raise ValueError('Cannot min-max normalize constant pooled time.')
        x['pseudotime']=(x.pseudotime-lo)/(hi-lo)
    elif lo < -1e-10 or hi > 1+1e-10:raise ValueError('Pseudotime must lie in [0,1]. Use normalize=True for a pooled min-max transform.')
    x['pseudotime']=x.pseudotime.clip(0,1)
    if not 0<bin_width<=1 or not np.isclose(1/bin_width,round(1/bin_width)):raise ValueError('bin_width must divide [0,1] exactly, e.g. 0.05.')
    vals=[x.loc[x.group==g,'pseudotime'].to_numpy() for g in groups];a,b=vals
    mw=stats.mannwhitneyu(b,a,alternative='two-sided',method='asymptotic',use_continuity=True)
    ks=stats.ks_2samp(a,b,alternative='two-sided',method='asymp')
    effect=float(2*mw.statistic/(len(a)*len(b))-1);stars=significance_stars(mw.pvalue)
    pooled=np.concatenate(vals)
    bw=float(bandwidth) if bandwidth is not None else max(float(np.std(pooled,ddof=1)*len(pooled)**(-.2)),.005)
    if not np.isfinite(bw) or bw<=0:raise ValueError('bandwidth must be finite and positive')
    prefix=Path(output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':11,'text.color':'#26364A','axes.labelcolor':'#26364A','axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'}):
        fig=plt.figure(figsize=(10.5,8.6));gs=fig.add_gridspec(2,1,height_ratios=[1,3.2],hspace=.22,left=.14,right=.91,top=.82,bottom=.18)
        top=fig.add_subplot(gs[0]);ax=fig.add_subplot(gs[1],sharex=top)
        bp=top.boxplot(vals,positions=[1,0],vert=False,widths=.48,patch_artist=True,showfliers=False,whis=1.5,medianprops={'color':'#172C3A','linewidth':2},whiskerprops={'color':'#526171'},capprops={'color':'#526171'})
        for patch,col in zip(bp['boxes'],colors):patch.set_facecolor(col);patch.set_alpha(.55);patch.set_edgecolor(col)
        top.set_yticks([1,0],[f'{g} (n={len(v)})' for g,v in zip(groups,vals)]);top.set_ylim(-.6,1.6);top.tick_params(axis='x',bottom=False,labelbottom=False);top.spines['bottom'].set_visible(False)
        # Axes-fraction X, data-coordinate Y: bracket stays outside the [0,1] axis.
        tr=top.get_yaxis_transform();top.plot([1.015,1.035,1.035,1.015],[0,0,1,1],transform=tr,clip_on=False,color='#26364A',lw=1.2)
        top.text(1.043,.5,stars,transform=tr,ha='left',va='center',fontsize=16,fontweight='bold',clip_on=False)
        grid=np.linspace(0,1,1001);edges=np.linspace(0,1,round(1/bin_width)+1);curves=[];bins=[]
        for g,v,col in zip(groups,vals,colors):
            density=np.zeros(len(grid))
            # Chunking keeps memory bounded for large cell counts.
            for start in range(0,len(v),1000):
                chunk=v[start:start+1000]
                for centers in [chunk,-chunk,2-chunk]:density+=np.exp(-.5*((grid[:,None]-centers[None,:])/bw)**2).sum(axis=1)
            density/=len(v)*bw*np.sqrt(2*np.pi);density/=np.trapezoid(density,grid)
            smooth=density*bin_width*100;counts=np.histogram(v,edges)[0];fractions=counts/len(v)*100
            ax.stairs(fractions,edges,color=col,linewidth=.9,alpha=.4);ax.fill_between(grid,smooth,color=col,alpha=.14);ax.plot(grid,smooth,color=col,lw=2.3,label=f'{g} · median {np.median(v):.3f}')
            curves.extend({'group':g,'pseudotime':t,'density':d,'percent_per_bin_width':s} for t,d,s in zip(grid,density,smooth))
            bins.extend({'group':g,'left':l,'right':r,'cells':int(n),'percent':f} for l,r,n,f in zip(edges[:-1],edges[1:],counts,fractions))
        ax.set(xlim=(0,1),ylim=(0,None),xlabel='Pseudotime on shared reference (0–1)',ylabel=f'Cells within group (% per {bin_width:g} pseudotime)');ax.set_xticks([0,.25,.5,.75,1]);ax.grid(axis='y',color='#E7ECF0',lw=.7);ax.set_axisbelow(True);ax.legend(frameon=False,loc='best')
        fig.text(.14,.95,title,fontsize=15,weight='bold')
        fig.text(.14,.905,f'KS: D = {ks.statistic:.3f}, P = {ks.pvalue:.2e}    |    Mann–Whitney U: P = {mw.pvalue:.2e}',fontsize=10)
        fig.text(.14,.872,f'Median shift ({groups[1]} − {groups[0]}) = {np.median(b)-np.median(a):+.3f}    |    Rank-biserial effect = {effect:+.3f}',fontsize=10,color='#56667A')
        fig.text(.14,.083,f'Thin steps: observed {bin_width:g}-width bins; curves: reflected KDE with shared bandwidth.',fontsize=8.5,color='#627185')
        fig.text(.14,.058,'Box: median / IQR, whiskers 1.5 × IQR. Stars: MWU raw P; * < .05, ** < .01, *** < .001.',fontsize=8.5,color='#627185')
        fig.text(.14,.033,'Cell-level exploratory comparison; not independent donor-level inference. Outliers omitted only in boxplot.',fontsize=8.5,color='#627185')
        paths={}
        for ext in ['jpg','png','svg']:
            dest=str(prefix)+'.'+ext;fig.savefig(dest,dpi=240,facecolor='white');paths[ext]=dest
        plt.close(fig)
    result={'group_col':group_col,'time_col':time_col,'groups':groups,'n_input':n_input,'excluded_invalid_or_unmapped':excluded_invalid,'excluded_unselected_groups':excluded_other,'normalize_pooled':normalize,'input_range':[lo,hi],'n':[len(a),len(b)],'medians':[float(np.median(v)) for v in vals],'MWU_p':float(mw.pvalue),'KS_D':float(ks.statistic),'KS_p':float(ks.pvalue),'rank_biserial_comparison_minus_reference':effect,'stars':stars,'star_test':'two-sided asymptotic MWU raw P, continuity and tie correction; no multiple-testing correction','bandwidth':bw,'bin_width':bin_width,'paths':paths}
    Path(str(prefix)+'_statistics.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    pd.DataFrame(curves).to_csv(str(prefix)+'_curves.csv',index=False);pd.DataFrame(bins).to_csv(str(prefix)+'_bins.csv',index=False);x.to_csv(str(prefix)+'_used_cells.csv',index=False)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input');ap.add_argument('--output',required=True);ap.add_argument('--group-col');ap.add_argument('--time-col');ap.add_argument('--groups',nargs=2);ap.add_argument('--normalize',action='store_true');ap.add_argument('--title',default='Pseudotime distributions');args=ap.parse_args()
    print(json.dumps(plot_pseudotime(args.input,args.output,args.group_col,args.time_col,args.groups,normalize=args.normalize,title=args.title),indent=2,ensure_ascii=False))
