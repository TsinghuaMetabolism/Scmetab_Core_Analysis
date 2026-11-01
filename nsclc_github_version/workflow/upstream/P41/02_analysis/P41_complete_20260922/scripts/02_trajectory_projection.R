suppressPackageStartupMessages({library(data.table);library(monocle);library(VGAM);library(Biobase);library(igraph)})
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));d<-readRDS(file.path(root,'audit/integrated.rds'));ref<-d$metadata$batch!='B3';cells<-as.data.frame(d$metadata[ref]);X<-d$after[ref,,drop=FALSE];ff<-as.data.frame(d$features)
N_ORDERING<-100L;NCENTER<-50L;SIGMA<-.1;SEED<-20260918L;KW_FDR<-.05;MIN_ABS_CLIFF<-.15;Z_CAP<-10
out<-function(x,n)fwrite(x,file.path(root,'panel_b/tables',n))
stage<-ifelse(cells$day=="D0","D0",ifelse(cells$day=="D3","D3","Endpoint"))
cliff<-function(a,b){r<-rank(c(a,b),ties.method="average");u<-sum(r[seq_along(a)])-length(a)*(length(a)+1)/2;1-2*u/(length(a)*length(b))}
evidence<-do.call(rbind,lapply(seq_len(ncol(X)),function(j){x<-X[,j];kw<-if(diff(range(x))>1e-12) kruskal.test(x,factor(stage)) else list(statistic=0,p.value=1);data.frame(feature_id=colnames(X)[j],KW_H=unname(kw$statistic),KW_p=kw$p.value,D0_endpoint_delta_pooled=cliff(x[stage=="D0"],x[stage=="Endpoint"]))}))
evidence$KW_FDR<-p.adjust(evidence$KW_p,"BH");evidence$KW_epsilon2<-pmax(0,(evidence$KW_H-2)/(nrow(X)-3))
evidence$eligible<-evidence$KW_FDR<KW_FDR & abs(evidence$D0_endpoint_delta_pooled)>=MIN_ABS_CLIFF
evidence<-evidence[order(!evidence$eligible,-evidence$KW_epsilon2,evidence$KW_FDR),];selected_ids<-head(evidence$feature_id[evidence$eligible],N_ORDERING)
# Preserve source matrix feature order, as in the original implementation.
selected_ids<-colnames(X)[colnames(X) %in% selected_ids]
Y<-X[,selected_ids,drop=FALSE];med<-apply(Y,2,median);mad0<-apply(abs(sweep(Y,2,med,"-")),2,median)
pop_sd<-apply(Y,2,function(v)sqrt(mean((v-mean(v))^2)))
scale0<-ifelse(mad0>1e-8,mad0,pop_sd+1e-8);Z<-sweep(sweep(Y,2,med,"-"),2,scale0,"/");Z[Z>Z_CAP]<-Z_CAP;Z[Z< -Z_CAP]<- -Z_CAP
evidence<-merge(evidence,ff[,c("feature_id","annotation_name","mz")],by="feature_id",sort=FALSE)
evidence$selected<-evidence$feature_id %in% selected_ids
out(evidence,'ordering_feature_screen.csv');sc<-data.frame(feature_id=selected_ids,reference_median=med,reference_scale_MAD=scale0);out(sc,'ordering_feature_scaling.csv')
fit_tree<-function(Z,ncenter=50,sigma=.1,seed=SEED){
 pd<-cells;rownames(pd)<-pd$cell_id;fd<-data.frame(gene_short_name=colnames(Z),row.names=colnames(Z))
 cds<-newCellDataSet(t(Z),phenoData=new('AnnotatedDataFrame',data=pd),featureData=new('AnnotatedDataFrame',data=fd),lowerDetectionLimit=-Inf,expressionFamily=VGAM::gaussianff())
 cds<-setOrderingFilter(cds,colnames(Z));set.seed(seed);cds<-reduceDimension(cds,max_components=2,reduction_method='DDRTree',norm_method='none',pseudo_expr=0,relative_expr=FALSE,auto_param_selection=FALSE,ncenter=ncenter,sigma=sigma,scaling=FALSE,verbose=FALSE)
 cds<-orderCells(cds);pd<-pData(cds);d0s<-table(pd$State[pd$day=='D0']);alls<-table(pd$State);frac<-d0s/alls[names(d0s)];cand<-names(d0s)[d0s>=max(5,ceiling(sum(pd$day=='D0')*.05))];rt<-as.numeric(cand[which.max(frac[cand])]);cds<-orderCells(cds,root_state=rt);cds
}
cds<-fit_tree(Z);pd<-pData(cds);xy<-t(reducedDimS(cds));cells$trajectory_1<-xy[cells$cell_id,1];cells$trajectory_2<-xy[cells$cell_id,2];cells$pseudotime<-pd[cells$cell_id,'Pseudotime'];cells$reference_Monocle_state<-as.character(pd[cells$cell_id,'State'])
g<-minSpanningTree(cds);ei<-as_edgelist(g,names=FALSE);kk<-reducedDimK(cds);nodes<-if(ncol(kk)==vcount(g))t(kk)else t(reducedDimS(cds));edges<-data.frame(x0=nodes[ei[,1],1],y0=nodes[ei[,1],2],x1=nodes[ei[,2],1],y1=nodes[ei[,2],2])
saveRDS(cds,file.path(root,'audit/Osi_Monocle2_model.rds'),compress=FALSE)
# Parameter checks are diagnostic only: never select by day correlation.
cfg<-data.frame(ncenter=c(50,50,30,70),sigma=c(.05,.2,.1,.1));checks<-list()
for(i in seq_len(nrow(cfg))){tt<-fit_tree(Z,cfg$ncenter[i],cfg$sigma[i]);pt<-pData(tt)[cells$cell_id,'Pseudotime'];checks[[i]]<-data.frame(cfg[i,],PT_spearman=cor(cells$pseudotime,pt,method='spearman'),n_states=length(unique(pData(tt)$State)),tree_branch_nodes=sum(degree(minSpanningTree(tt))>2))}
out(rbindlist(checks),'parameter_stability.csv');out(data.table(n_ordering=length(selected_ids),n_states=length(unique(cells$reference_Monocle_state)),branch_nodes=sum(degree(g)>2)),'reference_topology.csv')
all<-as.data.frame(d$metadata);all$trajectory_1<-NA_real_;all$trajectory_2<-NA_real_;all$pseudotime<-NA_real_;all$reference_Monocle_state<-NA_character_
all[ref,c('trajectory_1','trajectory_2','pseudotime','reference_Monocle_state')]<-cells[,c('trajectory_1','trajectory_2','pseudotime','reference_Monocle_state')]
XS<-sweep(sweep(d$after[,selected_ids,drop=FALSE],2,med,'-'),2,pmax(scale0,1e-6),'/');XR<-XS[ref,,drop=FALSE];XQ<-XS[!ref,,drop=FALSE];Y<-as.matrix(cells[,c('trajectory_1','trajectory_2','pseudotime')])
ridge_fit<-function(x,y){xm<-colMeans(x);ym<-colMeans(y);xc<-sweep(x,2,xm,'-');yc<-sweep(y,2,ym,'-');coef<-solve(crossprod(xc)+diag(10,ncol(x)),crossprod(xc,yc));list(coefficients=coef,intercept=ym-as.vector(xm%*%coef))}
predict_ridge<-function(m,x)sweep(x%*%m$coefficients,2,m$intercept,'+')
model<-ridge_fit(XR,Y);predraw<-predict_ridge(model,XQ);shift<-apply(predraw[all$day[!ref]=='D0',,drop=FALSE],2,median)-apply(Y[cells$day=='D0',,drop=FALSE],2,median);pred<-sweep(predraw,2,shift,'-');all[!ref,colnames(Y)]<-pred
r<-rank(all$pseudotime[ref],ties.method='first');all$progress_state<-NA_character_;all$progress_state[ref]<-as.character(cut(r,quantile(r,c(0,1/3,2/3,1)),include.lowest=TRUE,labels=c('Early','Middle','Late')));threshold<-quantile(cells$pseudotime,c(1/3,2/3));all$progress_state[!ref]<-as.character(cut(pred[,3],c(-Inf,threshold,Inf),labels=c('Early','Middle','Late')))
all$trajectory_coordinate_source<-ifelse(ref,'Osi_DDRTree_reference','fixed_reference_ridge_D0_calibrated');all$regimen<-all$treatment_series
# Stratified held-out cell CV, and reference-nearest-neighbor support in ordering-feature space.
set.seed(17);fold<-integer(sum(ref));for(s in unique(interaction(cells$batch,cells$day))){ix<-which(interaction(cells$batch,cells$day)==s);fold[ix]<-sample(rep(1:5,length.out=length(ix)))}
cv<-Y*NA;for(k in 1:5)cv[fold==k,]<-predict_ridge(ridge_fit(XR[fold!=k,,drop=FALSE],Y[fold!=k,,drop=FALSE]),XR[fold==k,,drop=FALSE])
fwrite(data.table(quantity=colnames(Y),RMSE=sqrt(colMeans((cv-Y)^2)),R2=1-colSums((cv-Y)^2)/colSums(sweep(Y,2,colMeans(Y),'-')^2)),file.path(root,'panel_g/tables/mapping_CV.csv'))
dd<-as.matrix(dist(XR));diag(dd)<-Inf;nr<-apply(dd,1,min);nq<-apply(XQ,1,function(x)min(sqrt(rowSums(sweep(XR,2,x,'-')^2))));cutoff<-quantile(nr,.95)
all$nearest_reference_distance<-NA_real_;all$nearest_reference_distance[ref]<-nr;all$nearest_reference_distance[!ref]<-nq;all$outside_reference_95pct<-all$nearest_reference_distance>cutoff
all$outside_reference_PT_range<-all$pseudotime<min(cells$pseudotime)|all$pseudotime>max(cells$pseudotime)
fwrite(data.table(parameter=colnames(Y),D0_translation=shift),file.path(root,'panel_g/tables/D0_calibration.csv'));fwrite(all[,c('cell_id','group','nearest_reference_distance','outside_reference_95pct','outside_reference_PT_range')],file.path(root,'panel_g/tables/off_reference_audit.csv'))
fwrite(all[ref,],file.path(root,'panel_b/inputs/Osi_cells_and_trajectory.csv'));fwrite(all,file.path(root,'panel_g/inputs/reference_and_projected_cells.csv'))
for(p in c('b','g'))fwrite(edges,file.path(root,paste0('panel_',p), 'inputs/fixed_reference_tree_edges.csv'))
d$metadata<-as.data.table(all);d$ordering<-sc;d$ridge<-model;d$ridge$D0_translation<-shift;saveRDS(d,file.path(root,'audit/analyzed.rds'),compress=FALSE)
print(table(all$group,all$progress_state));print(rbindlist(checks))
