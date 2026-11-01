src<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(src))))
suppressPackageStartupMessages(library(ggplot2))
a<-read.csv('inputs/selected_25_GAM_both_scenarios.csv')
cur<-read.csv('inputs/fitted_curves_95CI_both_scenarios.csv')
names_map<-read.csv('inputs/name_mapping.csv')
a$short_name<-names_map$short_name[match(a$feature_id,names_map$feature_id)]
stopifnot(all(cur$scenario=='exclude_D0'))
groups<-list(two=c('F00039','F02419'),seven=c('F07676','F08244','F07438','F08765','F08062','F07704','F08274'))
pal<-c('#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#56B4E9','#6A51A3')
wrap<-function(x,width=100){
  if(nchar(x)<=width)return(x)
  tok<-strsplit(gsub('-','- ',x,fixed=TRUE),' +')[[1]];out<-character();line<-''
  for(w in tok){sep<-if(nzchar(line)&&!grepl('-$',line))' ' else ''
    if(nchar(paste0(line,sep,w))>width&&nzchar(line)){out<-c(out,line);line<-w}else line<-paste0(line,sep,w)}
  paste(c(out,line),collapse='\n')
}
audits<-list()
for(g in names(groups)){
  ids<-groups[[g]];aa<-a[match(ids,a$feature_id),];stopifnot(!anyNA(aa$annotation))
  d<-cur[cur$feature_id%in%ids,];d$condition<-factor(d$condition,levels=c('Osi','Osi+TCD'));d$feature_id<-factor(d$feature_id,levels=ids)
  stopifnot(nrow(d)==length(ids)*500,all(table(d$feature_id,d$condition)==250))
  colors<-setNames(pal[seq_along(ids)],ids);labs<-setNames(vapply(aa$short_name,wrap,character(1)),ids)
  p<-ggplot(d,aes(t,fit,colour=feature_id,fill=feature_id,group=feature_id))+
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.14,colour=NA)+geom_line(linewidth=.85)+
    facet_wrap(~condition,nrow=1,scales='fixed',axes='all',axis.labels='all')+
    scale_colour_manual(values=colors,labels=labs)+scale_fill_manual(values=colors,labels=labs)+
    guides(colour=guide_legend(ncol=1),fill=guide_legend(ncol=1))+
    scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(0,.005)))+
    scale_y_continuous(expand=expansion(mult=c(.1,.15)))+
    labs(title=if(g=='two')'Choline and GPC' else 'Selected choline-related lipids',
      x='Pseudotime',y='Relative intensity',colour=NULL,fill=NULL)+
    theme_classic(base_size=14,base_family='Arial')+
    theme(plot.title=element_text(size=13,face='bold',hjust=.5,lineheight=1.05),
      axis.title=element_text(size=15,face='bold'),axis.title.x=element_text(margin=margin(t=10)),
      axis.title.y=element_text(margin=margin(r=8)),axis.text=element_text(colour='black',size=12),
      axis.line=element_line(linewidth=.65),axis.ticks=element_line(linewidth=.65),
      strip.background=element_blank(),strip.text=element_text(size=12,face='bold'),
      panel.spacing=grid::unit(1.4,'lines'),legend.position='bottom',legend.text=element_text(size=10),
      legend.key.width=grid::unit(.75,'cm'),legend.key.height=grid::unit(.45,'cm'),
      legend.margin=margin(t=10),plot.margin=margin(14,20,8,12))
  name<-if(g=='two')'Choline_and_GPC' else 'Seven_choline_related_lipids'
  height<-if(g=='two')5.8 else 7.1
  ggsave(paste0('figures/',name,'.pdf'),p,width=10.4,height=height,device=cairo_pdf,bg='white')
  ggsave(paste0('preview/',name,'.png'),p,width=10.4,height=height,dpi=180,bg='white')
  aa$group<-g;aa$color<-colors;write.csv(aa,paste0('tables/',g,'_metabolite_index.csv'),row.names=FALSE)
  for(id in ids)for(cond in levels(d$condition)){
    z<-d[d$feature_id==id & d$condition==cond,];z<-z[order(z$t),]
    audits[[length(audits)+1]]<-data.frame(group=g,feature_id=id,annotation=aa$annotation[match(id,ids)],condition=cond,
      t_start=z$t[1],t_end=tail(z$t,1),fit_start=z$fit[1],fit_end=tail(z$fit,1),
      fitted_endpoint_change=tail(z$fit,1)-z$fit[1],fitted_range=max(z$fit)-min(z$fit))
  }
}
write.csv(do.call(rbind,audits),'tables/descriptive_curve_changes.csv',row.names=FALSE)
