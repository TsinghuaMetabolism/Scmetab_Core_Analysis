POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.82;POINT_STROKE<-.18;POINT_BORDER<-'black'
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(source_file))))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
d<-fread('inputs/Osi_cells_and_trajectory.csv')
e<-fread('inputs/fixed_reference_tree_edges.csv')
s<-fread('inputs/ordering_feature_screen.csv')
stopifnot(nrow(d)==1068L,nrow(e)==49L,sum(s$selected)==159L,
          uniqueN(d$Monocle_state)==1L)
pal<-c(untreated_D0='#2166AC',Osi_D3='#2CB1BC',Osi_D8='#F6AD31',Osi_D14='#D1493F')
d[,group:=factor(group,levels=names(pal))]
# Rotate only the display because the fitted curve is long and narrow.
b<-ggplot(d,aes(trajectory_2,trajectory_1))+
  geom_segment(data=e,aes(x=y0,y=x0,xend=y1,yend=x1),inherit.aes=FALSE,
               linewidth=.55,colour='#647887')+
  coord_equal()+theme_classic(base_size=12)+
  theme(legend.position='bottom',plot.title.position='plot',
        plot.title=element_text(face='bold'))+
  labs(x='DDRTree 2 (rotated display)',y='DDRTree 1')
point<-function(a)geom_point(a,shape=POINT_SHAPE,size=POINT_SIZE,
  alpha=POINT_ALPHA,stroke=POINT_STROKE,colour=POINT_BORDER)
p<-(b+point(aes(fill=group))+scale_fill_manual(values=pal,name='Sampling day')+
      guides(fill=guide_legend(nrow=2))+
      labs(title='P39 Osimertinib trajectory',
           subtitle='159 D0/D14-significant features | ncenter = 50 | sigma = 0.1'))+
   (b+point(aes(fill=pseudotime))+scale_fill_viridis_c(option='C',name='Pseudotime')+
      labs(title='D0-rooted pseudotime',subtitle='One natural Monocle state'))
ggsave('figures/panel_b_Osi_reference_trajectory.pdf',p,width=14,height=6,
       device=cairo_pdf,bg='white')
ggsave('figures/panel_b_Osi_reference_trajectory.png',p,width=14,height=6,
       dpi=300,bg='white')
fwrite(d,'tables/plot_cells.csv');fwrite(e,'tables/plot_edges.csv')
