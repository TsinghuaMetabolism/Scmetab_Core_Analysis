POINT_SHAPE<-21;POINT_SIZE<-2;POINT_ALPHA<-.85;POINT_STROKE<-.2;POINT_BORDER<-'black'
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
read_input<-function(name)read.csv(file.path('inputs',name),check.names=FALSE,stringsAsFactors=FALSE)
read_matrix<-function(name){x<-read_input(name);ids<-x[[1]];m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-ids;storage.mode(m)<-'double';m}
suppressPackageStartupMessages(library(ggplot2))
d<-list(cells=read_input('reference_and_projected_cells.csv'),edges=read_input('fixed_reference_tree_edges.csv'));cells<-d$cells;e<-d$edges;poly<-read_input('reference_local_support_polygons.csv');q<-cells[cells$batch=='B3',];ref<-cells[cells$batch!='B3',]
pal<-c(untreated_D0='#2166AC',OsiTCD_D3='#7251A1',OsiTCD_D8='#A23B72');q$group<-factor(q$group,levels=names(pal))
theme_set(theme_classic(base_size=12)+theme(legend.position='bottom',plot.title.position='plot',plot.caption.position='plot',plot.title=element_text(face='bold')))
point<-function(data,a=NULL,...)geom_point(data=data,mapping=a,shape=POINT_SHAPE,size=POINT_SIZE,alpha=POINT_ALPHA,stroke=POINT_STROKE,colour=POINT_BORDER,...)
p<-ggplot(q,aes(trajectory_1,trajectory_2))+geom_polygon(data=poly,aes(x,y,group=polygon_id),inherit.aes=FALSE,fill='#D6E9ED',alpha=.18,colour=NA)+point(ref,fill='#D7DADD')+geom_segment(data=e,aes(x=x0,y=y0,xend=x1,yend=y1),inherit.aes=FALSE,linewidth=.5)+point(q,aes(fill=group))+scale_fill_manual(values=pal,name='Projected query')+coord_equal()+labs(x='Fixed DDRTree 1',y='Fixed DDRTree 2',title='Osi+TCD projected onto the fixed Osi reference',subtitle='Grey: fixed Osi cells; colour: constrained Ridge query\nLocal supports follow each P41 reference-tree arm')
savep<-function(p,name,w,h){ggsave(paste0(name,'.pdf'),p,width=w,height=h,device=cairo_pdf,bg='white');ggsave(paste0(name,'.png'),p,width=w,height=h,dpi=300,bg='white')}
xy_counts<-aggregate(rep(1L,nrow(q)),q[c('trajectory_1','trajectory_2')],sum);names(xy_counts)[3]<-'n';pile<-xy_counts[which.max(xy_counts$n),]
p<-p+labs(caption=sprintf('%d projected cells at %d distinct coordinates; %d coincide at one boundary vertex.',nrow(q),nrow(xy_counts),pile$n))
if(pile$n>1L)p<-p+annotate('segment',x=pile$trajectory_1-1.45,y=pile$trajectory_2-.28,xend=pile$trajectory_1-.14,yend=pile$trajectory_2-.025,colour='#404040',linewidth=.35,arrow=grid::arrow(length=grid::unit(.1,'inches')))+annotate('label',x=pile$trajectory_1-1.55,y=pile$trajectory_2-.35,label=sprintf('n = %d',pile$n),hjust=1,size=3.2,fill='white',label.size=.2)
savep(p,'figures/panel_g1_fixed_reference_projection',9,6.8)
cells$regimen<-ifelse(cells$batch=='B3','Osi+TCD','Osi');states<-c('Early','Middle','Late')
rows<-list();tests<-list()
for(day in c('D3','D8')){
 z<-cells[cells$day==day,];tab<-table(factor(z$regimen,levels=c('Osi','Osi+TCD')),factor(z$progress_state,levels=states))
 tab_test<-tab[,colSums(tab)>0,drop=FALSE];cs<-suppressWarnings(chisq.test(tab_test,correct=FALSE)); small<-any(cs$expected<5);pv<-if(small)fisher.test(tab_test)$p.value else cs$p.value;tests[[day]]<-data.frame(day=day,test=if(small)'Fisher exact' else 'Pearson chi-square',chi2=unname(cs$statistic),p_value=pv,cramers_v=sqrt(unname(cs$statistic)/sum(tab)),minimum_expected_count=min(cs$expected))
 v<-as.data.frame(tab);names(v)<-c('regimen','state','n');v$day<-day;v$percent<-100*v$n/ave(v$n,v$regimen,FUN=sum);rows[[day]]<-v
}
ct<-do.call(rbind,rows);ts<-do.call(rbind,tests);ts$BH_FDR<-p.adjust(ts$p_value,'BH');ct$state<-factor(ct$state,levels=states)
label<-setNames(sprintf('%s | P = %.3g; Cramer\'s V = %.3f',ts$day,ts$p_value,ts$cramers_v),ts$day)
p<-ggplot(ct,aes(regimen,percent,fill=state))+geom_col(width=.65,colour='white',linewidth=.35,position=position_stack(reverse=TRUE))+geom_text(aes(label=sprintf('%.1f%%\n(n=%d)',percent,n)),position=position_stack(vjust=.5,reverse=TRUE),size=3.6)+facet_wrap(~day,nrow=1,labeller=labeller(day=label))+scale_fill_manual(values=c(Early='#6EA6BC',Middle='#E4B354',Late='#C9494F'),name='Reference interval')+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=expansion(mult=c(0,.02)))+labs(x=NULL,y='Cells (%)',title='Matched-day redistribution along the Osi reference',subtitle='Reference-only graph-distance K=3 intervals; descriptive, not biological states')
savep(p,'figures/panel_g2_state_composition',11,5.5)
write.csv(cells,'tables/plot_cells.csv',row.names=FALSE);write.csv(e,'tables/plot_edges.csv',row.names=FALSE);write.csv(ct,'tables/composition.csv',row.names=FALSE);write.csv(ts,'tables/composition_tests.csv',row.names=FALSE)
