suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
a<-grep('^--file=',commandArgs(),value=TRUE)
setwd(dirname(dirname(normalizePath(sub('^--file=','',a[1])))))
d<-fread('inputs/Osi_cells_and_trajectory.csv')
e<-fread('inputs/fixed_reference_tree_edges.csv')
days<-c('D0','D3','D8','D14');stages<-c('Early','Middle','Late')
d[,day:=factor(day,levels=days)];d[,progress_state:=factor(progress_state,levels=stages)]
day_pal<-c(D0='#2166AC',D3='#2CB1BC',D8='#F6AD31',D14='#D1493F')
stage_pal<-c(Early='#3977B5',Middle='#E7AB3B',Late='#C55253')
base<-ggplot(d,aes(trajectory_2,trajectory_1))+
  geom_segment(data=e,aes(x=y0,y=x0,xend=y1,yend=x1),inherit.aes=FALSE,
               colour='#a9b5bf',linewidth=.4)+coord_equal()+theme_classic(base_size=11)+
  labs(x='DDRTree 2 (rotated display)',y='DDRTree 1')
p_batch<-base+geom_point(aes(fill=batch),shape=21,size=1.65,stroke=.15,alpha=.8)+
  scale_fill_manual(values=c(B1='#276FBF',B2='#E87524'))+
  labs(title='Batch distribution')+theme(legend.position='bottom')
p_pt<-ggplot(d,aes(day,pseudotime,fill=batch))+
  geom_violin(position=position_dodge(.8),scale='width',alpha=.55,trim=FALSE)+
  geom_boxplot(position=position_dodge(.8),width=.15,outlier.shape=NA)+
  scale_fill_manual(values=c(B1='#276FBF',B2='#E87524'))+theme_classic(base_size=11)+
  labs(title='Pseudotime by actual day',x='Actual day',y='Pseudotime',fill='Batch')+
  theme(legend.position='bottom')
p<-p_batch+p_pt+plot_annotation(title='P39 | Batch and true-time distribution on the new trajectory',
  subtitle='D3 is B1-only; D8 is B2-only')
ggsave('figures/panel_b_batch_and_time.pdf',p,width=13,height=5.8,device=cairo_pdf,bg='white')
ggsave('figures/panel_b_batch_and_time.png',p,width=13,height=5.8,dpi=300,bg='white')

tab<-d[,.(n=.N),by=.(day,progress_state)][order(day,progress_state)]
tab[,day_n:=sum(n),by=day];tab[,pct_within_day:=100*n/day_n]
fwrite(tab,'tables/stage_distribution_by_day.csv')
q<-ggplot(tab,aes(day,pct_within_day,fill=progress_state))+
  geom_col(width=.68,position=position_stack(reverse=TRUE))+
  geom_text(aes(label=sprintf('%d\n%.1f%%',n,pct_within_day)),
            position=position_stack(vjust=.5,reverse=TRUE),size=3.4,fontface='bold')+
  scale_fill_manual(values=stage_pal)+scale_y_continuous(labels=function(x)paste0(x,'%'))+
  theme_classic(base_size=12)+theme(legend.position='bottom')+
  labs(title='Early / Middle / Late composition by sampling day',
       subtitle='Global pseudotime-rank thirds; descriptive display stages',
       x='Actual day',y='Percent within day',fill='Stage')
ggsave('figures/panel_b_stage_composition.pdf',q,width=8,height=5,device=cairo_pdf,bg='white')
ggsave('figures/panel_b_stage_composition.png',q,width=8,height=5,dpi=300,bg='white')
fwrite(d[,.(cell_id,sample_id,group,batch,day,pseudotime,Monocle_state,progress_state)],
       'tables/trajectory_cell_summary.csv')
