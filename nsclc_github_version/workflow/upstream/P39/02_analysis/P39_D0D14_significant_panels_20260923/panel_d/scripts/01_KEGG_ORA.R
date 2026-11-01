# KEGG ORA from the exact annotated feature set displayed in panel C.
DISPLAY_NOMINAL_P<-.05
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE)
setwd(dirname(dirname(normalizePath(sub('^--file=','',a[1])))))
f<-fread('inputs/features.csv');s<-fread('inputs/all_features_GAM_significance.csv')
path<-unique(fread('inputs/TERM2GENE.csv'));pn<-fread('inputs/TERM2NAME.csv')
stopifnot(nrow(s)==3338L,identical(f$feature_id,s$feature_id))
q<-merge(f[!is.na(annotation)&nzchar(trimws(annotation))],
         s[association_BH_FDR_Osi<.05],by='feature_id',sort=FALSE)
fwrite(q,'tables/GAM_significant_annotated_query_features.csv')
ids<-function(x){
 if(is.na(x)||!nzchar(x))return(character())
 z<-unlist(strsplit(x,'[/;]'))
 unique(z[grepl('^C[0-9]{5}$',z)])
}
full<-sort(unique(path$KEGGID))
det<-intersect(unique(unlist(lapply(f$KEGG_IDs,ids))),full)
ora<-function(universe,label){
 query<-intersect(unique(unlist(lapply(q$KEGG_IDs,ids))),universe)
 N<-length(universe);n<-length(query)
 r<-rbindlist(lapply(unique(path$PathwayID),function(id){
   bg<-intersect(unique(path[PathwayID==id,KEGGID]),universe);M<-length(bg)
   if(M<2L)return(NULL)
   hits<-intersect(bg,query);k<-length(hits)
   hf<-q[vapply(KEGG_IDs,function(x)any(ids(x)%in%hits),logical(1))]
   data.table(ID=id,Description=pn$PathwayName[match(id,pn$PathwayID)],Count=k,
     query_n=n,pathway_n=M,universe_N=N,
     pvalue=if(n)phyper(k-1,M,N-M,n,lower.tail=FALSE)else 1,
     FoldEnrich=if(n)(k/n)/(M/N)else NA_real_,KEGG_IDs=paste(hits,collapse='/'),
     Feature_IDs=paste(hf$feature_id,collapse=';'))
 }))
 r[,p.adjust:=NA_real_];r[Count>0,p.adjust:=p.adjust(pvalue,'BH')]
 r[,p.adjust_all_eligible:=p.adjust(pvalue,'BH')]
 fwrite(r,file.path('tables',paste0('ORA_',label,'.csv')))
 fwrite(data.table(KEGGID=query),file.path('tables',paste0('query_',label,'.csv')))
 fwrite(data.table(KEGGID=universe),file.path('tables',paste0('background_',label,'.csv')))
 r
}
r<-ora(full,'all_pathway');ora(det,'detected')
z<-r[Count>0&pvalue<DISPLAY_NOMINAL_P]
setorder(z,p.adjust,-FoldEnrich,Description)
z[,top_to_bottom_order:=seq_len(.N)]
fwrite(z,'tables/plot_data.csv')
z[,Description:=factor(Description,levels=rev(Description))]
p<-ggplot(z,aes(-log10(p.adjust),Description))+
  geom_segment(aes(x=0,xend=-log10(p.adjust),yend=Description),
               colour='#CCD5DF',linewidth=.7)+
  geom_point(aes(size=Count,fill=FoldEnrich),shape=21,colour='black',
             stroke=.2,alpha=.87)+
  geom_vline(xintercept=-log10(.05),linetype=2,colour='#B2182B')+
  scale_size_area(max_size=9)+
  scale_fill_gradientn(colours=c('#6AADC8','#F6D57E','#D1493F','#A50026'))+
  theme_classic(base_size=11)+
  theme(axis.title.y=element_blank(),plot.title.position='plot',
        plot.caption.position='plot',plot.title=element_text(face='bold'))+
  labs(x=expression(-log[10]('BH-FDR')),fill='Fold enrichment',
       size='Unique KEGG\ncompounds',
       title='P39 | KEGG enrichment of pseudotime-associated metabolites',
       subtitle='Full pathway-compound background | nominal P < 0.05 displayed',
       caption='Top: smaller BH-FDR; ties: larger fold enrichment. Dashed line: BH-FDR = 0.05.')
h<-max(5,.30*nrow(z)+2)
ggsave('figures/panel_d_KEGG_enrichment.pdf',p,width=12,height=h,device=cairo_pdf,bg='white')
ggsave('figures/panel_d_KEGG_enrichment.png',p,width=12,height=h,dpi=300,bg='white')
cat('Query features:',nrow(q),'unique compounds:',nrow(fread('tables/query_all_pathway.csv')),
    'displayed pathways:',nrow(z),'\n')
