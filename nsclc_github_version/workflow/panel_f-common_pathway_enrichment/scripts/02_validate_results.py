#!/usr/bin/env python3
from pathlib import Path
import csv,hashlib,math,json,itertools
root=Path(__file__).resolve().parents[1]
def read(p):
 with (root/p).open(newline='') as f:return list(csv.DictReader(f))
def ids(s):return {i for i in s.replace(';','/').split('/') if len(i)==6 and i.startswith('C') and i[1:].isdigit()}
def bh(p):
 order=sorted(range(len(p)),key=p.__getitem__);out=[None]*len(p);cur=1
 for rank in range(len(p),0,-1):
  i=order[rank-1];cur=min(cur,p[i]*len(p)/rank);out[i]=cur
 return out
files=read('inputs/source_manifest.csv')
for r in files:assert hashlib.sha256((root/r['copied_path']).read_bytes()).hexdigest()==r['sha256']
features=read('inputs/joint_trajectory/shared_features_ordered.csv')
mat={}
for d in ['P39','P41','P59']:
 mat[d]={r['shared_id']:[float(r[f'T{i:03d}']) for i in range(1,161)] for r in read(f'inputs/joint_trajectory/{d}_shared_rowZ.csv')}
expected_two=set();expected_three=set()
for f in features:
 pairs=[]
 for a,b in itertools.combinations(['P39','P41','P59'],2):
  corr=sum(x*y for x,y in zip(mat[a][f['shared_id']],mat[b][f['shared_id']]))/159
  assert abs(corr-float(f[f'cor_{a}_{b}']))<1e-12
  if float(f[a+'_q'])<.05 and float(f[b+'_q'])<.05 and f[a+'_shape']==f[b+'_shape']!='Complex' and corr>=.5:pairs.append((a,b))
 all_corr=[sum(x*y for x,y in zip(mat[a][f['shared_id']],mat[b][f['shared_id']]))/159 for a,b in itertools.combinations(['P39','P41','P59'],2)]
 shapes=[f[d+'_shape'] for d in ['P39','P41','P59']]
 n_sig=sum(float(f[d+'_q'])<.05 for d in ['P39','P41','P59'])
 if n_sig>=2 and len(set(shapes))==1 and shapes[0]!='Complex' and min(all_corr)>=.5:expected_two.add(f['shared_id'])
 if len(pairs)==3:expected_three.add(f['shared_id'])
assert (len(expected_three),len(expected_two))==(20,38)
heatmap_ids={r['shared_id'] for r in read('inputs/joint_trajectory/common_all_three_trends_at_least_two_significant.csv')}
assert expected_two==heatmap_ids
sets={}
for r in read('inputs/P59_panel_d_reference/inputs/TERM2GENE.csv'):sets.setdefault(r['PathwayID'],set()).add(r['KEGGID'])
universe=set().union(*sets.values());N=len(universe);assert N==2901
maxerr=0;summary={}
for name,expected in [('at_least_two_donors',expected_two)]:
 q=read(f'results/tables/{name}_metabolite_list.csv');assert {r['shared_id'] for r in q}==expected
 query=set().union(*(ids(r['KEGG_IDs']) for r in q))&universe;n=len(query)
 rows=read(f'results/tables/{name}_KEGG_all_eligible.csv');assert len(rows)==86
 ps=[]
 for r in rows:
  pathway=sets[r['ID']];M=len(pathway);k=len(pathway&query)
  p=sum(math.comb(M,j)*math.comb(N-M,n-j) for j in range(k,min(n,M)+1) if 0<=n-j<=N-M)/math.comb(N,n)
  assert int(r['Count'])==k and ids(r['KEGG_IDs'])==pathway&query
  assert abs(float(r['FoldEnrich'])-(k/n)/(M/N))<1e-10
  assert abs(float(r['pvalue'])-p)<1e-12
  maxerr=max(maxerr,abs(float(r['pvalue'])-p));ps.append(p)
 hit=[i for i,r in enumerate(rows) if int(r['Count'])>0];adjust=bh([ps[i] for i in hit]);alladjust=bh(ps)
 for i,qval in zip(hit,adjust):assert abs(float(rows[i]['p.adjust'])-qval)<1e-12
 for r,qval in zip(rows,alladjust):assert abs(float(r['p.adjust_all_eligible'])-qval)<1e-12
 summary[name]=dict(features=len(q),unique_compounds=n,hit_BH_significant=sum(x<.05 for x in adjust),all_eligible_BH_significant=sum(x<.05 for x in alladjust))
report=dict(passed=True,verified_input_files=len(files),selection_scope='91 features shared by all three donors',background_compounds=N,ORA_max_absolute_pvalue_difference=maxerr,lists=summary)
(root/'results/validation/independent_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
