# Stratified 80/20 validation: no held-out intensities enter affine estimation or PCA fitting.
SEED<-20260922L;TRAIN_FRACTION<-.8;COMBAT_TOL<-1e-4
suppressPackageStartupMessages({library(data.table);library(cluster)})
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));setwd(root)
d<-readRDS('data/integration_candidates.rds');m<-copy(d$meta);X<-d$before;stopifnot(nrow(m)==1068,!anyNA(m$day),all(m$day!='Unknown'),identical(rownames(X),m$cell_id))
out<-file.path(root,'heldout_validation');dir.create(out,showWarnings=FALSE)
# Import estimator definitions only; do not run or overwrite the main integration.
expr<-parse('scripts/01_integrate.R');for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c('combat_affine','apply_affine'))eval(e)
anchor_days<-c('D0','D14')
set.seed(SEED);m[,partition:='not_in_split']
for(b in unique(m$batch))for(day0 in anchor_days){
  ix<-which(m$batch==b&m$day==day0)
  m[ix,partition:='heldout']
  m[sample(ix,floor(TRAIN_FRACTION*length(ix))),partition:='train']
}
train<-m$partition=='train';test<-m$partition=='heldout';not_split<-m$partition=='not_in_split'
stopifnot(all(m$day[train|test]%in%anchor_days),all(m$day[not_split]%in%c('D3','D8')),!any(m$cell_id[train]%in%m$cell_id[test]))
fwrite(m[,.(cell_id,sample_id,batch,day,partition)],file.path(out,'cell_partitions.csv'));counts<-m[,.N,by=.(batch,day,partition)];fwrite(counts,file.path(out,'partition_counts.csv'))
tr0<-combat_affine(X,m$batch,train&m$day=='D0');base<-apply_affine(X,m$batch,tr0);tr14<-combat_affine(base,m$batch,train&m$day=='D14');second<-apply_affine(base,m$batch,tr14);w<-as.numeric(sub('D','',m$day))/14;dual<-base+(second-base)*w
# Evaluation coordinate system learned from D0/D14 training cells only.
pc<-prcomp(X[train,,drop=FALSE],center=TRUE,scale.=FALSE,rank.=50);proj<-function(x)sweep(x,2,pc$center,'-')%*%pc$rotation;pcs<-lapply(list(Before=X,D0_only=base,D0_D14=dual),proj)
knn<-function(z,k=15L){dd<-as.matrix(dist(z));diag(dd)<-Inf;t(apply(dd,1,function(v)order(v)[seq_len(min(k,nrow(z)-1))]))}
# Five balanced resamplings ONLY within fixed held-out cells: not five independent validation splits.
set.seed(SEED+1L);metrics<-list()
for(day0 in c('D0','D14'))for(rep in 1:5){ix<-which(test&m$day==day0);n<-min(table(m$batch[ix]));take<-unlist(lapply(c('B1','B2'),function(b)sample(ix[m$batch[ix]==b],n)));bb<-m$batch[take]
for(nm in names(pcs)){z<-pcs[[nm]][take,,drop=FALSE];nn<-knn(z);p<-rowMeans(matrix(bb[nn],nrow(nn))=='B1');metrics[[length(metrics)+1]]<-data.table(anchor=day0,repeat_id=rep,method=nm,n_per_batch=n,local_mixing=mean(1/(p^2+(1-p)^2)),mean_abs_batch_ASW=mean(abs(silhouette(as.integer(factor(bb)),dist(z))[,3])),centroid_gap=sqrt(sum((colMeans(z[bb=='B1',])-colMeans(z[bb=='B2',]))^2)))}}
metrics<-rbindlist(metrics);fwrite(metrics,file.path(out,'heldout_balanced_metrics.csv'));sm<-metrics[,lapply(.SD,mean),by=.(anchor,method,n_per_batch),.SDcols=c('local_mixing','mean_abs_batch_ASW','centroid_gap')];fwrite(sm,file.path(out,'heldout_metric_means.csv'));fwrite(metrics[,lapply(.SD,sd),by=.(anchor,method,n_per_batch),.SDcols=c('local_mixing','mean_abs_batch_ASW','centroid_gap')],file.path(out,'heldout_metric_subsampling_sd.csv'))
bio<-rbindlist(lapply(c('D0_only','D0_D14'),function(nm)rbindlist(lapply(c('B1','B2'),function(b){ix<-which(test&m$batch==b);z<-pcs$Before[ix,];v<-pcs[[nm]][ix,];n0<-knn(z);n1<-knn(v);jac<-sapply(seq_along(ix),function(i)length(intersect(n0[i,],n1[i,]))/length(union(n0[i,],n1[i,])));data.table(method=nm,batch=b,n_heldout=length(ix),distance_spearman=cor(as.vector(dist(z)),as.vector(dist(v)),method='spearman'),knn_jaccard=mean(jac))}))));fwrite(bio,file.path(out,'heldout_biological_structure.csv'))
delta<-rbindlist(lapply(c('D0_only','D0_D14'),function(nm)rbindlist(lapply(c('B1','B2'),function(b){i0<-test&m$batch==b&m$day=='D0';i14<-test&m$batch==b&m$day=='D14';a<-colMeans(X[i14,])-colMeans(X[i0,]);v<-if(nm=='D0_only')base else dual;z<-colMeans(v[i14,])-colMeans(v[i0,]);data.table(method=nm,batch=b,delta_spearman=cor(a,z,method='spearman'),direction_agreement=mean(sign(a)==sign(z)),delta_norm_ratio=sqrt(sum(z^2)/sum(a^2)))}))));fwrite(delta,file.path(out,'heldout_D0_D14_conservation.csv'))
params<-rbindlist(lapply(c('D0','D14'),function(st){tt<-if(st=='D0')tr0 else tr14;rbindlist(lapply(names(tt),function(b)data.table(stage=st,batch=b,feature_id=colnames(X),slope=tt[[b]]$slope,offset=tt[[b]]$offset,iterations=tt[[b]]$iterations)))}));fwrite(params,file.path(out,'training_only_affine_parameters.csv'));saveRDS(list(seed=SEED,training_ids=m$cell_id[train],heldout_ids=m$cell_id[test],not_in_split_ids=m$cell_id[not_split],PCA=pc,D0=tr0,D14=tr14),file.path(out,'training_fitted_transforms.rds'),compress=FALSE)
# Independent re-estimation on physically subsetted training matrices audits the fitting masks.
r0<-combat_affine(X[train,,drop=FALSE],m$batch[train],m$day[train]=='D0');r14<-combat_affine(base[train,,drop=FALSE],m$batch[train],m$day[train]=='D14');err<-max(abs(unlist(tr0)-unlist(r0)),abs(unlist(tr14)-unlist(r14)))
anchor_check<-m[day%in%anchor_days,.(ok=sum(partition=='train')==floor(.8*.N)&&sum(partition=='heldout')==.N-floor(.8*.N)),by=.(batch,day)]
checks<-data.table(
  check=c('anchor train and heldout disjoint','all anchor cells assigned train or heldout','D3/D8 excluded from split','training fraction floor(0.8*n) in each anchor stratum','physically training-only refit equals masked fit','PCA uses anchor training cells only'),
  passed=c(!any(m$cell_id[train]%in%m$cell_id[test]),sum(train)+sum(test)==sum(m$day%in%anchor_days),all(m$partition[m$day%in%c('D3','D8')]=='not_in_split'),all(anchor_check$ok),err<1e-12,nrow(pc$x)==sum(train)),
  detail=c(paste(sum(train),'anchor train;',sum(test),'anchor heldout'),paste(sum(m$day%in%anchor_days),'anchor cells'),paste(sum(not_split),'D3/D8 cells'),'',paste('max parameter error',err),paste(nrow(pc$x),'anchor training rows')))
fwrite(checks,file.path(out,'leakage_checks.csv'));stopifnot(all(checks$passed))
print(counts);print(sm);print(bio);print(delta);capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
# Uncompressed outputs allow an independent reviewer to recompute metrics without original raw files.
for(nm in c('Before','D0_only','D0_D14')){mat<-switch(nm,Before=X,D0_only=base,D0_D14=dual);fwrite(data.table(cell_id=m$cell_id[test],as.data.table(mat[test,,drop=FALSE])),file.path(out,paste0('heldout_',nm,'_intensity.csv')));fwrite(data.table(cell_id=m$cell_id[test],as.data.table(pcs[[nm]][test,,drop=FALSE])),file.path(out,paste0('heldout_',nm,'_PC50.csv')))}
fwrite(m[test],file.path(out,'heldout_cell_metadata.csv'))
suppressPackageStartupMessages(library(ggplot2));long<-melt(sm,id.vars=c('anchor','method','n_per_batch'));long[,method:=factor(method,c('Before','D0_only','D0_D14'))];long[,variable:=factor(variable,levels=c('local_mixing','mean_abs_batch_ASW','centroid_gap'),labels=c('Local mixing (higher is better)','Absolute batch ASW (lower is better)','Centroid gap (lower is better)'))]
p<-ggplot(long,aes(method,value,group=anchor,colour=anchor))+geom_line()+geom_point(aes(fill=anchor),shape=21,size=2,alpha=.85,stroke=.2,colour='black')+geom_text(aes(label=sprintf('%.3f',value),vjust=ifelse(anchor=='D0',-1,1.6)),size=3,show.legend=FALSE)+facet_wrap(~variable,scales='free_y',ncol=3)+scale_colour_manual(values=c(D0='#2166AC',D14='#D1493F'))+scale_fill_manual(values=c(D0='#2166AC',D14='#D1493F'))+scale_y_continuous(expand=expansion(mult=c(.14,.2)))+theme_classic(base_size=11)+theme(legend.position='bottom',axis.text.x=element_text(angle=20,hjust=1))+labs(title='P39 | Anchor-training correction evaluated on held-out anchors',subtitle=paste(sum(train),'D0/D14 training cells;',sum(test),'D0/D14 held-out cells;',sum(not_split),'D3/D8 cells not in split'),x=NULL,y=NULL,caption='Stratified 80/20 split only within D0 and D14. PCA and correction parameters use anchor training cells only.')
ggsave(file.path(out,'heldout_evaluation.pdf'),p,width=13,height=5.5,device=cairo_pdf);ggsave(file.path(out,'heldout_evaluation.png'),p,width=13,height=5.5,dpi=180)
