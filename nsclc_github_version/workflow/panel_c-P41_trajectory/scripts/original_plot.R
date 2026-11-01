POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.85;POINT_STROKE<-.2;POINT_BORDER<-'black'
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
read_input<-function(name)read.csv(file.path('inputs',name),check.names=FALSE,stringsAsFactors=FALSE)
read_matrix<-function(name){x<-read_input(name);ids<-x[[1]];m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-ids;storage.mode(m)<-'double';m}
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
x<-list(cells=read_input('Osi_cells_and_trajectory.csv'),edges=read_input('fixed_reference_tree_edges.csv'));d<-x$cells;e<-x$edges
pal<-c(untreated_D0='#2166AC',Osi_D3='#2CB1BC',Osi_D8='#F6AD31',Osi_D14='#D1493F');d$group<-factor(d$group,levels=names(pal))
mid<-c(mean(range(d$trajectory_1)),mean(range(d$trajectory_2)));span<-max(diff(range(d$trajectory_1)),diff(range(d$trajectory_2)))*.55
b<-ggplot(d,aes(trajectory_1,trajectory_2))+geom_segment(data=e,aes(x=x0,y=y0,xend=x1,yend=y1),inherit.aes=FALSE,linewidth=.5)+coord_equal(xlim=mid[1]+c(-span,span),ylim=mid[2]+c(-span,span))+theme_classic(base_size=12)+theme(legend.position='bottom',plot.title.position='plot',plot.title=element_text(face='bold'))+labs(x='DDRTree 1',y='DDRTree 2')
point<-function(a)geom_point(a,shape=POINT_SHAPE,size=POINT_SIZE,alpha=POINT_ALPHA,stroke=POINT_STROKE,colour=POINT_BORDER)
p<-(b+point(aes(fill=group))+scale_fill_manual(values=pal,name='Treatment / time')+guides(fill=guide_legend(nrow=2))+labs(title='Osi monotherapy',subtitle='100 features | ncenter = 50 | sigma = 0.1'))+(b+point(aes(fill=pseudotime))+scale_fill_viridis_c(option='C',name='Pseudotime')+labs(title='Reference pseudotime',subtitle='Fixed Osi reference; native states retained in CSV'))
ggsave('figures/panel_b_Osi_reference_trajectory.pdf',p,width=13,height=5.6,device=cairo_pdf,bg='white');ggsave('figures/panel_b_Osi_reference_trajectory.png',p,width=13,height=5.6,dpi=300,bg='white')
write.csv(d,'tables/plot_cells.csv',row.names=FALSE);write.csv(e,'tables/plot_edges.csv',row.names=FALSE)
