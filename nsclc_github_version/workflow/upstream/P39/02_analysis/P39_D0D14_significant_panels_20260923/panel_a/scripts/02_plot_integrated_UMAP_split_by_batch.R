# Show B1 and B2 in the same jointly integrated UMAP coordinate system.
POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.85;POINT_STROKE<-.2;POINT_BORDER<-'black'
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(source_file))))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table)})
d<-fread('inputs/cells_and_UMAP_coordinates.csv')
stopifnot(nrow(d)==1068L,!anyDuplicated(d$cell_id),all(d$batch%in%c('B1','B2')),
          !anyNA(d$cell_type),all(nzchar(d$cell_type)))
d[,`:=`(x=UMAP_after_1,y=UMAP_after_2)]
limits_x<-range(d$x);limits_y<-range(d$y);tic_limits<-range(log10(d$TIC))
pal<-c(untreated_D0='#2166AC',Osi_D3='#2CB1BC',Osi_D8='#F6AD31',Osi_D14='#D1493F')
plots<-list()
for(v in c('group','TIC'))for(batch0 in c('B1','B2')){
 z<-d[batch==batch0]
 z[,display:=if(v=='TIC')log10(TIC) else factor(group,levels=names(pal))]
 p<-ggplot(z,aes(x,y,fill=display))+
   geom_point(shape=POINT_SHAPE,size=POINT_SIZE,alpha=POINT_ALPHA,
              stroke=POINT_STROKE,colour=POINT_BORDER)+
   coord_equal(xlim=limits_x,ylim=limits_y,expand=TRUE)+theme_classic(base_size=11)+
   theme(plot.title=element_text(face='bold'),plot.title.position='plot')+
   labs(x='Integrated UMAP 1',y='Integrated UMAP 2',
        title=paste0(batch0,' | ',v,' | n = ',nrow(z)),
        fill=if(v=='TIC')'log10 TIC' else 'Time')
 if(v=='group')p<-p+scale_fill_manual(values=pal,drop=TRUE)
 if(v=='TIC')p<-p+scale_fill_viridis_c(option='C',limits=tic_limits)
 plots[[length(plots)+1L]]<-p
}
p<-wrap_plots(plots,ncol=2)+
  plot_annotation(title='P39 | Batch-specific views of the integrated UMAP',
    subtitle='Same corrected coordinates, axis limits and TIC color scale; no re-embedding',
    caption='D3 occurs only in B1; D8 occurs only in B2.')
ggsave('figures/panel_a_integrated_UMAP_split_by_batch.pdf',p,width=14,height=9,
       device=cairo_pdf,bg='white')
ggsave('figures/panel_a_integrated_UMAP_split_by_batch.png',p,width=14,height=9,
       dpi=300,bg='white')
for(batch0 in c('B1','B2'))fwrite(d[batch==batch0],
  file.path('tables',paste0('integrated_UMAP_',batch0,'_cells.csv')))
check<-d[,.(n_cells=.N,UMAP1_min=min(x),UMAP1_max=max(x),
            UMAP2_min=min(y),UMAP2_max=max(y)),by=batch]
fwrite(check,'tables/split_UMAP_validation.csv')
