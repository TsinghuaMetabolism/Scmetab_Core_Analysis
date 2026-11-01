src<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(src))))
suppressPackageStartupMessages(library(ggplot2))
st<-read.csv('tables/curve_difference_statistics.csv')
cur<-read.csv('tables/fitted_curves_95CI_both_scenarios.csv');cur<-cur[cur$scenario=='include_D0',]
cur$condition<-factor(cur$condition,levels=c('Osi','Osi+TCD'))
stopifnot(nrow(st)==25,!anyDuplicated(st$pdf_filename))
pal<-c('Osi'='#2384BB','Osi+TCD'='#E58B12')
wrap_title<-function(s) {
  s<-gsub('/',' / ',s,fixed=TRUE)
  if(nchar(s)<=42)return(s)
  # Split at chemical hyphens while preserving the complete annotation.
  tokens<-strsplit(gsub('-','- ',s,fixed=TRUE),' +')[[1]]
  lines<-character();line<-''
  for(w in tokens){sep<-if(nchar(line)>0 && !grepl('-$',line))' ' else ''
    if(nchar(paste0(line,sep,w))>43 && nzchar(line)){lines<-c(lines,line);line<-w}else line<-paste0(line,sep,w)}
  paste(c(lines,line),collapse='\n')
}
for(scenario in c('include_D0','exclude_D0')) {
st<-read.csv('tables/selected_25_GAM_both_scenarios.csv');st<-st[st$scenario==scenario,]
cur<-read.csv('tables/fitted_curves_95CI_both_scenarios.csv');cur<-cur[cur$scenario==scenario,]
cur$condition<-factor(cur$condition,levels=c('Osi','Osi+TCD'))
plots<-list()
for(i in seq_len(nrow(st))) {
  d<-cur[cur$feature_id==st$feature_id[i],]
  qt<-formatC(st$overall_BH_FDR[i],format='e',digits=2)
  assoc_label<-function(g,q)paste0(g,'  ',if(q<.001)'***' else if(q<.01)'**' else if(q<.05)'*' else 'ns',
                                '  (trend q = ',formatC(q,format='e',digits=2),')')
  leg<-c(assoc_label('Osi',st$association_BH_FDR_Osi[i]),assoc_label('Osi+TCD',st$association_BH_FDR_TCDOsi[i]))
  p<-ggplot(d,aes(t,fit,colour=condition,fill=condition))+
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.14,colour=NA)+geom_line(linewidth=.85)+
    scale_colour_manual(values=pal,labels=leg)+scale_fill_manual(values=pal,labels=leg)+
    guides(colour=guide_legend(ncol=1),fill=guide_legend(ncol=1))+
    scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(0,.005)))+
    scale_y_continuous(expand=expansion(mult=c(.1,.15)))+
    labs(title=wrap_title(st$annotation[i]),subtitle=paste0('Curve difference: q = ',qt),caption=paste0(if(scenario=='include_D0')'D0 included' else 'D0 excluded (sensitivity)',' | GAM | BH: 1,928 features'),
         x='Pseudotime',y='Relative intensity',colour=NULL,fill=NULL)+
    theme_classic(base_size=14,base_family='Arial')+
    theme(plot.title=element_text(size=13,face='bold',hjust=.5,lineheight=1.05),
      plot.subtitle=element_text(size=10,hjust=.5,margin=margin(b=10)),
      axis.title=element_text(size=15,face='bold'),axis.title.x=element_text(margin=margin(t=10)),
      axis.title.y=element_text(margin=margin(r=8)),axis.text=element_text(colour='black',size=12),
      axis.line=element_line(linewidth=.65),axis.ticks=element_line(linewidth=.65),
      plot.caption=element_text(size=8,hjust=.5,colour='#555555'),legend.position='bottom',legend.text=element_text(size=10,face='bold'),
      legend.key.width=grid::unit(.65,'cm'),plot.margin=margin(14,20,8,12))
  plots[[i]]<-p
  if(scenario=='include_D0') {
  ggsave(file.path('figures',st$pdf_filename[i]),p,width=5.2,height=5.3,device=cairo_pdf,bg='white')
  ggsave(file.path('preview',paste0(st$feature_id[i],'.png')),p,width=5.2,height=5.3,dpi=180,bg='white')
  }
}
cairo_pdf(if(scenario=='include_D0')'All_25_metabolites.pdf' else 'Sensitivity_exclude_D0_25_metabolites.pdf',width=5.2,height=5.3,onefile=TRUE);for(p in plots)print(p);dev.off()
}
capture.output(sessionInfo(),file='tables/R_sessionInfo.txt')
