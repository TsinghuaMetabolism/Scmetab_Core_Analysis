# Adjustable cell point style; all inputs are CSV files in ../inputs.
POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.85;POINT_STROKE<-.2;POINT_BORDER<-'black'
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
read_input<-function(name)read.csv(file.path('inputs',name),check.names=FALSE,stringsAsFactors=FALSE)
read_matrix<-function(name){x<-read_input(name);ids<-x[[1]];m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-ids;storage.mode(m)<-'double';m}
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
d<-read_input('cells_and_UMAP_coordinates.csv')
pal<-c(untreated_D0='#2166AC',Osi_D3='#2CB1BC',Osi_D8='#F6AD31',Osi_D14='#D1493F',OsiTCD_D3='#7251A1',OsiTCD_D8='#A23B72')
d$group<-factor(d$group,levels=names(pal))
plots<-list();tab<-list()
for(v in c('group','batch','TIC'))for(s in c('before','after')){
 x<-d;x$x<-x[[paste0('UMAP_',s,'_1')]];x$y<-x[[paste0('UMAP_',s,'_2')]];x$display<-if(v=='TIC')log10(x$TIC) else x[[v]]
 p<-ggplot(x,aes(x,y,fill=display))+geom_point(shape=POINT_SHAPE,size=POINT_SIZE,alpha=POINT_ALPHA,stroke=POINT_STROKE,colour=POINT_BORDER)+coord_equal()+theme_classic(base_size=10)+theme(plot.title=element_text(face='bold'),plot.title.position='plot')+labs(x='UMAP 1',y='UMAP 2',title=paste(if(s=='before')'Before correction' else 'D0 + Osi D14 correction',v,sep=' | '),fill=if(v=='TIC')'log10 TIC' else v)
 if(v=='group')p<-p+scale_fill_manual(values=pal,drop=FALSE)
 if(v=='batch')p<-p+scale_fill_manual(values=c(B1='#276FBF',B2='#E87524',B3='#7251A1'))
 if(v=='TIC')p<-p+scale_fill_viridis_c(option='C')
 plots[[length(plots)+1]]<-p
}
p<-wrap_plots(plots,ncol=2)+plot_annotation(title='P41 | Batch integration',subtitle='Both embeddings use TIC normalization + log1p; frozen final coordinates')
ggsave('figures/panel_a_UMAP_before_after.pdf',p,width=14,height=14,device=cairo_pdf,bg='white');ggsave('figures/panel_a_UMAP_before_after.png',p,width=14,height=14,dpi=300,bg='white')
write.csv(d,'tables/plot_data.csv',row.names=FALSE)
