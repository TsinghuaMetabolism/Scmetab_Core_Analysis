# Primary sort: BH-FDR ascending (stronger at top); ties: Fold enrichment descending.
DISPLAY_NOMINAL_P<-.05
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
read_input<-function(name)read.csv(file.path('inputs',name),check.names=FALSE,stringsAsFactors=FALSE)
read_matrix<-function(name){x<-read_input(name);ids<-x[[1]];m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-ids;storage.mode(m)<-'double';m}
suppressPackageStartupMessages(library(ggplot2))
# Derive the query from the exact same Gaussian-GAM Osi association test as panel c/h.
s<-read_input('all_features_GAM_significance.csv')
a<-read_input('annotated_feature_KEGG_pathway_mapping.csv')
a$association_BH_FDR_Osi<-s$association_BH_FDR_Osi[match(a$feature_id,s$feature_id)]
a<-a[!is.na(a$association_BH_FDR_Osi)&a$association_BH_FDR_Osi<.05,]
path<-unique(read_input('TERM2GENE.csv'));pn<-read_input('TERM2NAME.csv')
split_ids<-function(x){if(is.na(x)||!nzchar(x))return(character());z<-unlist(strsplit(x,'[/;]'));unique(z[grepl('^C[0-9]{5}$',z)])}
universe<-sort(unique(path$KEGGID));query<-intersect(unique(unlist(lapply(a$KEGG_IDs,split_ids))),universe)
N<-length(universe);n<-length(query)
sets<-split(path$KEGGID,path$PathwayID)
d<-do.call(rbind,lapply(names(sets),function(id){
  bg<-unique(sets[[id]]);M<-length(bg);hits<-intersect(bg,query);k<-length(hits)
  if(M<2)return(NULL)
  hf<-a[vapply(a$KEGG_IDs,function(x)any(split_ids(x)%in%hits),logical(1)),]
  data.frame(ID=id,Description=pn$PathwayName[match(id,pn$PathwayID)],Count=k,
    GeneRatio=paste0(k,'/',n),BgRatio=paste0(M,'/',N),
    pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),FoldEnrich=(k/n)/(M/N),
    KEGG_IDs=paste(hits,collapse='/'),Feature_IDs=paste(hf$feature_id,collapse=';'))
}))
hit<-d$Count>0;d$p.adjust<-NA_real_;d$p.adjust[hit]<-p.adjust(d$pvalue[hit],'BH')
d$p.adjust_all_eligible<-p.adjust(d$pvalue,'BH')
write.csv(d,'tables/KEGG_enrichment_all_eligible_pathways.csv',row.names=FALSE)
write.csv(a,'tables/GAM_significant_annotated_query_features.csv',row.names=FALSE)
write.csv(data.frame(KEGGID=query),'tables/query_unique_KEGG_compounds.csv',row.names=FALSE)
write.csv(data.frame(KEGGID=universe),'tables/background_all_pathway_compounds.csv',row.names=FALSE)
d<-d[d$Count>0&d$pvalue<DISPLAY_NOMINAL_P,]
d<-d[order(d$p.adjust,-d$FoldEnrich,d$Description),];d$top_to_bottom_order<-seq_len(nrow(d))
write.csv(d,'tables/plot_data.csv',row.names=FALSE)
d$Description<-factor(d$Description,levels=rev(d$Description))
p<-ggplot(d,aes(-log10(p.adjust),Description))+geom_segment(aes(x=0,xend=-log10(p.adjust),yend=Description),colour='#CCD5DF',linewidth=.7)+geom_point(aes(size=Count,fill=FoldEnrich),shape=21,colour='black',stroke=.2,alpha=.85)+geom_vline(xintercept=-log10(.05),linetype=2,colour='#B2182B')+scale_size_area(max_size=9,breaks=sort(unique(d$Count)))+scale_fill_gradientn(colours=c('#6AADC8','#F6D57E','#D1493F','#A50026'))+theme_classic(base_size=11)+theme(axis.title.y=element_blank(),plot.title.position='plot',plot.caption.position='plot',plot.title=element_text(face='bold'))+labs(x=expression(-log[10]('BH-FDR')),fill='Fold enrichment',size='Unique KEGG\ncompounds',title='P41 | KEGG enrichment of pseudotime-associated metabolites',subtitle='Full pathway background | displayed pathways: nominal P < 0.05',caption='Top: smaller BH-FDR; ties: larger fold enrichment. Dashed line: BH-FDR = 0.05.')
ggsave('figures/panel_d_KEGG_enrichment.pdf',p,width=12,height=6,device=cairo_pdf,bg='white');ggsave('figures/panel_d_KEGG_enrichment.png',p,width=12,height=6,dpi=300,bg='white')
