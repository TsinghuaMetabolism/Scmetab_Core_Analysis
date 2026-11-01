# D0/D14-only ordering-feature selection requested by the user. Require an
# endpoint difference in both batches with concordant direction; D3/D8 never
# participate in feature screening, but remain in the fitted trajectory.
suppressPackageStartupMessages({
  library(data.table); library(monocle); library(VGAM); library(igraph)
  library(ggplot2); library(patchwork)
})
a <- grep('^--file=', commandArgs(), value=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=', '', a[1]))))
setwd(root)
out <- 'trajectory_review/anchor_consistent_sensitivity'
dir.create(out, recursive=TRUE, showWarnings=FALSE)

d <- readRDS('data/final_analysis.rds')
m <- as.data.frame(d$meta)
X <- d$X
stopifnot(identical(d$choice,'D0_D14'), identical(rownames(X),m$cell_id))
cliff <- function(a,b) {
  r <- rank(c(a,b))
  u <- sum(r[seq_along(a)]) - length(a)*(length(a)+1)/2
  1 - 2*u/(length(a)*length(b))
}

bb <- rbindlist(lapply(c('B1','B2'), function(batch0) {
  i <- m$batch==batch0 & m$day=='D0'
  j <- m$batch==batch0 & m$day=='D14'
  delta <- vapply(seq_len(ncol(X)), function(k) cliff(X[i,k],X[j,k]), numeric(1))
  p <- vapply(seq_len(ncol(X)), function(k)
    wilcox.test(X[i,k],X[j,k],exact=FALSE)$p.value, numeric(1))
  data.table(feature_id=colnames(X),batch=batch0,cliff=delta,
             endpoint_FDR=p.adjust(p,'BH'))
}))
feat <- dcast(bb,feature_id~batch,value.var=c('cliff','endpoint_FDR'))
feat[, same_direction := sign(cliff_B1)==sign(cliff_B2) & sign(cliff_B1)!=0]
feat[, eligible := same_direction & abs(cliff_B1)>=.15 & abs(cliff_B2)>=.15 &
       endpoint_FDR_B1<.05 & endpoint_FDR_B2<.05]
feat[, weakest_batch_effect := pmin(abs(cliff_B1),abs(cliff_B2))]
feat[, worst_batch_FDR := pmax(endpoint_FDR_B1,endpoint_FDR_B2)]
setorder(feat,-eligible,-weakest_batch_effect,worst_batch_FDR,feature_id)
chosen <- feat[eligible==TRUE,feature_id]
stopifnot(length(chosen)>=3L)
ids <- colnames(X)[colnames(X)%in%chosen]
feat[, selected := feature_id %in% ids]
fwrite(feat,file.path(out,'feature_selection.csv'))

Y <- X[,ids,drop=FALSE]
med <- apply(Y,2,median)
mad0 <- apply(abs(sweep(Y,2,med,'-')),2,median)
sd0 <- apply(Y,2,function(z) sqrt(mean((z-mean(z))^2)))
sc <- ifelse(mad0>1e-8,mad0,sd0+1e-8)
Z <- sweep(sweep(Y,2,med,'-'),2,sc,'/')
Z[Z>10] <- 10; Z[Z< -10] <- -10
pd <- m; rownames(pd) <- pd$cell_id
fd <- data.frame(gene_short_name=ids,row.names=ids)
cds <- newCellDataSet(t(Z),
  phenoData=new('AnnotatedDataFrame',data=pd),
  featureData=new('AnnotatedDataFrame',data=fd),
  lowerDetectionLimit=-Inf,expressionFamily=gaussianff())
cds <- setOrderingFilter(cds,ids)
set.seed(20260918L)
cds <- reduceDimension(cds,max_components=2,reduction_method='DDRTree',
  norm_method='none',pseudo_expr=0,relative_expr=FALSE,
  auto_param_selection=FALSE,ncenter=50,sigma=.1,scaling=FALSE,verbose=FALSE)
cds <- orderCells(cds)
root_pd <- pData(cds)
tb <- table(root_pd$State[root_pd$day=='D0'])
tt <- table(root_pd$State)
fraction <- tb/tt[names(tb)]
ca <- names(tb)[tb>=max(5,ceiling(sum(root_pd$day=='D0')*.05))]
root_state <- as.numeric(ca[which.max(fraction[ca])])
cds <- orderCells(cds,root_state=root_state)
pd_out <- pData(cds)
xy <- t(reducedDimS(cds))
z <- data.table(cell_id=m$cell_id,day=m$day,batch=m$batch,
  trajectory_1=xy[m$cell_id,1],trajectory_2=xy[m$cell_id,2],
  pseudotime=pd_out[m$cell_id,'Pseudotime'],
  state=as.character(pd_out[m$cell_id,'State']))
fwrite(z,file.path(out,'cell_trajectory.csv'))
g <- minSpanningTree(cds)
ee <- as_edgelist(g,names=FALSE)
kk <- reducedDimK(cds)
nodes <- if(ncol(kk)==vcount(g)) t(kk) else t(reducedDimS(cds))
edges <- data.table(x0=nodes[ee[,1],1],y0=nodes[ee[,1],2],
                    x1=nodes[ee[,2],1],y1=nodes[ee[,2],2])
fwrite(edges,file.path(out,'trajectory_edges.csv'))
saveRDS(cds,file.path(out,'monocle2_model.rds'),compress=FALSE)

trend <- fread('trajectory_review/ordering_time_trend_audit.csv')
main_sel <- fread('tables/ordering_feature_screen.csv')[selected==TRUE,feature_id]
main_m <- as.data.table(d$meta)
quant <- function(v,p) as.numeric(quantile(v,p,names=FALSE))
pt_summary <- z[, .(n=.N,median_pseudotime=median(pseudotime),
                    q25=quant(pseudotime,.25),q75=quant(pseudotime,.75)),
                by=.(batch,day)][order(match(day,c('D0','D3','D8','D14')),batch)]
fwrite(pt_summary,file.path(out,'pseudotime_by_sample.csv'))
metrics <- data.table(
  model=c('previous_four_day_KW_100','D0_D14_significant_159'),
  n_features=c(100L,length(ids)),
  n_pooled_monotonic=c(sum(trend[feature_id%in%main_sel,pooled_monotonic]),
                       sum(trend[feature_id%in%ids,pooled_monotonic])),
  n_same_endpoint_direction=c(sum(trend[feature_id%in%main_sel,endpoint_same_direction]),
                              sum(trend[feature_id%in%ids,endpoint_same_direction])),
  n_within_batch_monotonic_both=c(sum(trend[feature_id%in%main_sel,within_batch_monotonic_both]),
                                  sum(trend[feature_id%in%ids,within_batch_monotonic_both])),
  day_pseudotime_spearman=c(
    cor(as.numeric(sub('D','',main_m$day)),main_m$pseudotime,method='spearman'),
    cor(as.numeric(sub('D','',z$day)),z$pseudotime,method='spearman')),
  B1_day_pseudotime_spearman=c(
    cor(as.numeric(sub('D','',main_m$day[main_m$batch=='B1'])),
        main_m$pseudotime[main_m$batch=='B1'],method='spearman'),
    cor(as.numeric(sub('D','',z$day[z$batch=='B1'])),
        z$pseudotime[z$batch=='B1'],method='spearman')),
  B2_day_pseudotime_spearman=c(
    cor(as.numeric(sub('D','',main_m$day[main_m$batch=='B2'])),
        main_m$pseudotime[main_m$batch=='B2'],method='spearman'),
    cor(as.numeric(sub('D','',z$day[z$batch=='B2'])),
        z$pseudotime[z$batch=='B2'],method='spearman'))
)
fwrite(metrics,file.path(out,'comparison.csv'))

pal <- c(D0='#2166AC',D3='#2CB1BC',D8='#F6AD31',D14='#D1493F')
z[,day:=factor(day,levels=names(pal))]
plot1 <- ggplot(z,aes(trajectory_2,trajectory_1,color=day))+
  geom_segment(data=edges,aes(x=y0,y=x0,xend=y1,yend=x1),
               inherit.aes=FALSE,color='#a8b4bd',linewidth=.35,alpha=.55)+
  geom_point(size=1.3,alpha=.73)+coord_equal()+
  scale_color_manual(values=pal)+theme_classic(base_size=11)+
  theme(legend.position='bottom')+
  labs(title='D0/D14-significant ordering features',
       subtitle=paste0(length(ids),' features significant in both batches; D3/D8 excluded from screening'),
       x='DDRTree 2',y='DDRTree 1',color='Day')
plot2 <- ggplot(z,aes(day,pseudotime,fill=batch))+
  geom_violin(position=position_dodge(width=.8),scale='width',alpha=.55,trim=FALSE)+
  geom_boxplot(position=position_dodge(width=.8),width=.16,outlier.shape=NA,alpha=.75)+
  scale_fill_manual(values=c(B1='#276FBF',B2='#D97524'))+
  theme_classic(base_size=11)+theme(legend.position='bottom')+
  labs(title='Pseudotime by actual day',x='Actual day',y='Pseudotime',fill='Batch')
plot <- plot1 / plot2
ggsave(file.path(out,'alternative_trajectory.png'),plot,
       width=11,height=10,dpi=180,bg='white')
ggsave(file.path(out,'alternative_trajectory.pdf'),plot,
       width=11,height=10,device=cairo_pdf,bg='white')
ggsave(file.path(out,'trajectory_by_day.png'),plot1,
       width=11,height=5,dpi=180,bg='white')
ggsave(file.path(out,'trajectory_by_day.pdf'),plot1,
       width=11,height=5,device=cairo_pdf,bg='white')
ggsave(file.path(out,'pseudotime_by_day.png'),plot2,
       width=9,height=5,dpi=180,bg='white')
ggsave(file.path(out,'pseudotime_by_day.pdf'),plot2,
       width=9,height=5,device=cairo_pdf,bg='white')

# Display the actual intermediate-day behavior of only the endpoint-selected
# features. The six columns remain batch-specific, avoiding a false pooled
# D3/D8 sequence when those two days occur in different batches.
sample_cols <- c('B1_D0','B1_D3','B1_D14','B2_D0','B2_D8','B2_D14')
sel_medians <- copy(trend[feature_id%in%ids,
  c('feature_id',sample_cols),with=FALSE])
sel_medians <- merge(sel_medians,feat[,.(feature_id,cliff_B1,cliff_B2)],
  by='feature_id',sort=FALSE)
fwrite(sel_medians,file.path(out,'selected_feature_sample_medians.csv'))
ord <- feat[eligible==TRUE][order(-sign(cliff_B1),-weakest_batch_effect),feature_id]
long <- melt(sel_medians,id.vars=c('feature_id','cliff_B1','cliff_B2'),
  measure.vars=sample_cols,variable.name='sample',value.name='median')
long[,batch:=sub('_.*','',as.character(sample))]
long[,z:= {s<-sd(median); if(s>0) (median-mean(median))/s else rep(0,.N)},
  by=.(feature_id,batch)]
long[,direction:=ifelse(cliff_B1>0,'D14 higher','D14 lower')]
long[,feature_id:=factor(feature_id,levels=rev(ord))]
long[,sample:=factor(sample,levels=sample_cols)]
heat <- ggplot(long,aes(sample,feature_id,fill=pmax(-2,pmin(2,z))))+
  geom_tile()+facet_grid(direction~.,scales='free_y',space='free_y')+
  geom_vline(xintercept=3.5,color='#304454',linewidth=.45)+
  scale_fill_gradient2(low='#2166AC',mid='#F7F7F7',high='#B2182B',
                       midpoint=0,name='Row Z')+
  theme_minimal(base_size=11)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle=35,hjust=1),panel.grid=element_blank())+
  labs(title='D0/D14-significant features across observed samples',
       subtitle='Within-batch row Z of corrected medians; D3 is B1-only, D8 is B2-only',
       x='Observed batch × day',y='159 selected features')
ggsave(file.path(out,'selected_feature_sample_medians.png'),heat,
       width=8,height=10,dpi=180,bg='white')
ggsave(file.path(out,'selected_feature_sample_medians.pdf'),heat,
       width=8,height=10,device=cairo_pdf,bg='white')

cat('Eligible:',sum(feat$eligible),'selected:',length(ids),
    'overlap with main:',length(intersect(ids,main_sel)),
    'states:',uniqueN(z$state),'root:',root_state,'\n')
print(metrics)
