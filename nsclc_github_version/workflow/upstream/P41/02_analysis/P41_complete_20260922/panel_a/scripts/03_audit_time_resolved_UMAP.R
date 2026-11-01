# Independent supplementary assessment from local panel CSV.
POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.85;POINT_STROKE<-.2;POINT_BORDER<-'black'
K_NEIGHBORS<-15L
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
d<-fread('inputs/cells_and_UMAP_coordinates.csv')
stopifnot(!anyDuplicated(d$cell_id),all(d$day%in%c('D0','D3','D8','D14')))
d[,cohort:=ifelse(batch=='B3','Osi+TCD','Osi')]
d[,display_group:=ifelse(day=='D0',paste0('Untreated D0 (',batch,')'),ifelse(batch=='B3',paste0('Osi+TCD ',day),paste0('Osi ',day,' (',batch,')')))]
counts<-d[,.(n=.N),by=.(day,batch,cohort,group,display_group)][order(match(day,c('D0','D3','D8','D14')),batch)]
fwrite(counts,'tables/time_resolved_cell_counts.csv')
long<-rbindlist(lapply(c('before','after'),function(stage){data.table(d[,.(cell_id,batch,day,cohort,group,display_group)],stage=stage,x=d[[paste0('UMAP_',stage,'_1')]],y=d[[paste0('UMAP_',stage,'_2')]])}))
long[,stage:=factor(stage,levels=c('before','after'),labels=c('Before correction','After D0 + Osi D14 correction'))]
long[,day:=factor(day,levels=c('D0','D3','D8','D14'))]
fwrite(long,'tables/time_resolved_UMAP_points.csv')
# Neighborhood diagnostics only compare same biological condition and day:
# untreated D0 across B1/B2/B3, and Osi D14 across B1/B2.
metrics<-list();percell<-list()
for(day0 in c('D0','D14'))for(stage0 in c('before','after')){
 z<-d[day==day0];xy<-as.matrix(z[,c(paste0('UMAP_',stage0,'_1'),paste0('UMAP_',stage0,'_2')),with=FALSE]);dd<-as.matrix(dist(xy));diag(dd)<-Inf
 near<-t(apply(dd,1,function(v)order(v)[1:K_NEIGHBORS]));bn<-matrix(z$batch[near],nrow=nrow(z),ncol=K_NEIGHBORS)
 lev<-sort(unique(z$batch));pr<-sapply(lev,function(b)rowMeans(bn==b));same<-rowMeans(bn==z$batch)
 score<-1/rowSums(pr^2)
 percell[[length(percell)+1]]<-data.table(cell_id=z$cell_id,day=day0,batch=z$batch,stage=stage0,same_batch_neighbor_fraction=same,local_inverse_simpson=score)
 metrics[[length(metrics)+1]]<-data.table(day=day0,stage=stage0,n_cells=nrow(z),n_batches=length(lev),mean_same_batch_neighbors=mean(same),mean_inverse_simpson=mean(score),median_inverse_simpson=median(score),fraction_all_15_same_batch=mean(same==1))
}
fwrite(rbindlist(percell),'tables/time_matched_UMAP_neighbor_mixing_per_cell.csv');fwrite(rbindlist(metrics),'tables/time_matched_UMAP_neighbor_mixing_summary.csv')
# D3/D8 compare different treatment conditions; overlap here is descriptive only.
contrast<-list()
for(day0 in c('D3','D8'))for(stage0 in c('before','after')){
 z<-d[day==day0];xy<-as.matrix(z[,c(paste0('UMAP_',stage0,'_1'),paste0('UMAP_',stage0,'_2')),with=FALSE]);dd<-as.matrix(dist(xy));diag(dd)<-Inf
 near<-t(apply(dd,1,function(v)order(v)[1:K_NEIGHBORS]));bn<-matrix(z$batch[near],nrow(z),K_NEIGHBORS)
 z[,other_batch_neighbor_fraction:=rowMeans(bn!=z$batch)]
 contrast[[length(contrast)+1]]<-z[,.(n=.N,mean_other_batch_neighbor_fraction=mean(other_batch_neighbor_fraction)),by=batch][,`:=`(day=day0,stage=stage0)]
}
fwrite(rbindlist(contrast),'tables/D3_D8_cross_treatment_neighbor_overlap.csv')
# Centroids and spread are within-embedding descriptions. Coordinate axes are independently fitted.
cent<-long[,.(n=.N,centroid_x=mean(x),centroid_y=mean(y),sd_x=sd(x),sd_y=sd(y),median_x=median(x),median_y=median(y)),by=.(stage,day,batch,cohort,group,display_group)]
fwrite(cent,'tables/time_resolved_UMAP_centroids.csv')
# Grey background conveys overall geometry; colored points are the current true-time slice.
background<-unique(long[,.(stage,x,y,cell_id)])
cols<-c(B1='#276FBF',B2='#E87524',B3='#7251A1')
xranges<-long[,.(xmin=min(x),xmax=max(x),ymin=min(y),ymax=max(y)),by=stage]
panel<-function(stage0,day0){z<-long[stage==stage0&day==day0];bg<-background[stage==stage0];rr<-xranges[stage==stage0]
 ggplot()+geom_point(data=bg,aes(x,y),shape=POINT_SHAPE,size=.7,alpha=.13,stroke=.08,colour='#B8C4CF',fill='#B8C4CF')+
 geom_point(data=z,aes(x,y,fill=batch),shape=POINT_SHAPE,size=POINT_SIZE,alpha=POINT_ALPHA,stroke=POINT_STROKE,colour=POINT_BORDER)+
 scale_fill_manual(values=cols,drop=FALSE)+coord_equal(xlim=c(rr$xmin,rr$xmax),ylim=c(rr$ymin,rr$ymax),expand=TRUE)+
 theme_classic(base_size=11)+theme(legend.position='bottom',plot.title=element_text(face='bold',size=12),axis.text=element_text(size=8))+
 labs(title=paste0(day0,' | ',as.character(stage0),' | n=',nrow(z)),x='UMAP 1',y='UMAP 2',fill='Batch')
}
suppressPackageStartupMessages(library(patchwork))
plots<-unlist(lapply(c('D0','D3','D8','D14'),function(day0)lapply(levels(long$stage),function(stage0)panel(stage0,day0))),recursive=FALSE)
p<-wrap_plots(plots,ncol=2)+plot_annotation(title='P41 | True-time cell distributions before and after batch correction',subtitle='Grey: all QC cells in that embedding. Colored: indicated day. D3/D8 include different treatments across batches.',caption='Each column shares limits across days; before/after UMAPs were fitted separately, so axis positions are not matched cell movements.')
ggsave('figures/panel_a_time_resolved_UMAP.pdf',p,width=14,height=24,device=cairo_pdf,bg='white')
ggsave('figures/panel_a_time_resolved_UMAP.png',p,width=14,height=24,dpi=220,bg='white')
print(rbindlist(metrics))
