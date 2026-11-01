# Subset the SAME jointly integrated coordinates. No new UMAP fit.
POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.85;POINT_STROKE<-.2;POINT_BORDER<-'black'
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
d<-read.csv('inputs/cells_and_UMAP_coordinates.csv',check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(d$cell_id),all(d$batch%in%c('B1','B2','B3')))
# Untreated D0 is assigned by experiment of origin, NOT by its common group label.
d$display_cohort<-ifelse(d$batch%in%c('B1','B2'),'Osi','Osi+TCD')
d$x<-d$UMAP_after_1;d$y<-d$UMAP_after_2
limits_x<-range(d$x);limits_y<-range(d$y);tic_limits<-range(log10(d$TIC))
pal<-c(untreated_D0='#2166AC',Osi_D3='#2CB1BC',Osi_D8='#F6AD31',Osi_D14='#D1493F',OsiTCD_D3='#7251A1',OsiTCD_D8='#A23B72')
batch_pal<-c(B1='#276FBF',B2='#E87524',B3='#7251A1')
plots<-list()
for(v in c('group','batch','TIC'))for(cohort in c('Osi','Osi+TCD')){
 z<-d[d$display_cohort==cohort,]
 z$display<-if(v=='TIC')log10(z$TIC) else if(v=='group')factor(z$group,levels=names(pal)) else z$batch
 p<-ggplot(z,aes(x,y,fill=display))+geom_point(shape=POINT_SHAPE,size=POINT_SIZE,alpha=POINT_ALPHA,stroke=POINT_STROKE,colour=POINT_BORDER)+coord_equal(xlim=limits_x,ylim=limits_y,expand=TRUE)+theme_classic(base_size=11)+theme(plot.title=element_text(face='bold'),plot.title.position='plot',legend.position='right')+labs(x='Integrated UMAP 1',y='Integrated UMAP 2',title=paste0(cohort,' | ',v,' | n = ',nrow(z)),fill=if(v=='TIC')'log10 TIC' else v)
 if(v=='group')p<-p+scale_fill_manual(values=pal,drop=TRUE)
 if(v=='batch')p<-p+scale_fill_manual(values=batch_pal,drop=TRUE)
 if(v=='TIC')p<-p+scale_fill_viridis_c(option='C',limits=tic_limits)
 plots[[length(plots)+1]]<-p
}
p<-wrap_plots(plots,ncol=2)+plot_annotation(title='P41 | Treatment-specific views of the integrated UMAP',subtitle='Same coordinates, axis limits and TIC color scale; no re-embedding',caption='Osi: B1 + B2, including their untreated D0 cells. Osi+TCD: B3, including its untreated D0 cells.')
ggsave('figures/panel_a_integrated_UMAP_split_by_treatment.pdf',p,width=14,height=13,device=cairo_pdf,bg='white')
ggsave('figures/panel_a_integrated_UMAP_split_by_treatment.png',p,width=14,height=13,dpi=300,bg='white')
cols<-c('sample_id','cell_id','group','batch','day','display_cohort','TIC','UMAP_after_1','UMAP_after_2')
write.csv(d[d$display_cohort=='Osi',cols],'tables/integrated_UMAP_Osi_cells.csv',row.names=FALSE)
write.csv(d[d$display_cohort=='Osi+TCD',cols],'tables/integrated_UMAP_OsiTCD_cells.csv',row.names=FALSE)
out<-rbind(read.csv('tables/integrated_UMAP_Osi_cells.csv'),read.csv('tables/integrated_UMAP_OsiTCD_cells.csv'))
idx<-match(out$cell_id,d$cell_id)
err<-max(abs(as.matrix(out[,c('UMAP_after_1','UMAP_after_2')])-as.matrix(d[idx,c('UMAP_after_1','UMAP_after_2')])))
stopifnot(nrow(out)==nrow(d),!anyDuplicated(out$cell_id),setequal(out$cell_id,d$cell_id),err<1e-12)
write.csv(data.frame(cohort=c('Osi','Osi+TCD'),n_cells=c(sum(d$display_cohort=='Osi'),sum(d$display_cohort=='Osi+TCD')),max_coordinate_difference=err),'tables/split_UMAP_validation.csv',row.names=FALSE)
