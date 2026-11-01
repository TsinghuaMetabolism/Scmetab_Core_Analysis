#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));setwd(root)
for(d in c('results/tables','results/figures','results/validation'))dir.create(d,recursive=TRUE,showWarnings=FALSE)
path<-unique(fread('inputs/P59_panel_d_reference/inputs/TERM2GENE.csv'))
pn<-fread('inputs/P59_panel_d_reference/inputs/TERM2NAME.csv')
p59ann<-fread('inputs/P59_panel_d_reference/inputs/annotated_feature_KEGG_pathway_mapping.csv',na.strings=c('','NA'))
map<-fread('inputs/joint_trajectory/shared_features_ordered.csv')
DONORS<-c('P39','P41','P59');FDR<-.05;COR<-.5
split_ids<-function(x){if(is.na(x)||!nzchar(x))return(character());z<-unlist(strsplit(x,'[/;]'));sort(unique(z[grepl('^C[0-9]{5}$',z)]))}
universe<-sort(unique(path$KEGGID));sets<-split(path$KEGGID,path$PathwayID)
fwrite(data.table(KEGGID=universe),'results/tables/background_all_pathway_compounds.csv')
# All list selection is restricted to the 91 matched features from panel_e.
pairs<-combn(DONORS,2,simplify=FALSE)
pa<-rbindlist(lapply(pairs,function(p){
 d1<-p[1];d2<-p[2];co<-map[[paste0('cor_',d1,'_',d2)]]
 q1<-map[[paste0(d1,'_q')]];q2<-map[[paste0(d2,'_q')]]
 s1<-map[[paste0(d1,'_shape')]];s2<-map[[paste0(d2,'_shape')]]
 data.table(shared_id=map$shared_id,annotation=map$annotation,pair=paste(p,collapse='+'),donor1=d1,donor2=d2,q1=q1,q2=q2,shape1=s1,shape2=s2,curve_correlation=co,both_significant=q1<FDR&q2<FDR,same_shape=s1==s2&s1!='Complex',qualifies=q1<FDR&q2<FDR&s1==s2&s1!='Complex'&co>=COR)
}));fwrite(pa,'results/tables/pairwise_selection_audit.csv')
qualified<-pa[qualifies==TRUE,.(qualifying_pairs=paste(pair,collapse=';'),n_qualifying_pairs=.N,qualifying_shapes=paste(sort(unique(shape1)),collapse=';')),by=shared_id]
map<-merge(map,qualified,by='shared_id',all.x=TRUE,sort=FALSE);map[is.na(n_qualifying_pairs),n_qualifying_pairs:=0L]
map[,selected_three:=concordant & all_three_significant]
map[,selected_at_least_two:=concordant & n_significant_donors>=2]
stopifnot(sum(map$selected_three)==20,sum(map$selected_at_least_two)==38,all(map[selected_three==TRUE]$selected_at_least_two))
heatmap_common<-fread('inputs/joint_trajectory/common_all_three_trends_at_least_two_significant.csv')
stopifnot(setequal(map[selected_at_least_two==TRUE]$shared_id,heatmap_common$shared_id))
fwrite(map[n_qualifying_pairs>0 & !selected_at_least_two],'results/tables/excluded_previous_pair_supported_only_features.csv')
map[,feature_id:=shared_id]
map[,KEGG_IDs:=p59ann$KEGG_IDs[match(P59_feature_id,p59ann$feature_id)]]
# Audit donor annotations; absent mappings are retained, conflicting nonempty IDs stop the run.
anns<-lapply(DONORS,function(d)fread(paste0('inputs/annotations/',d,'.csv'),na.strings=c('','NA')));names(anns)<-DONORS
ma<-rbindlist(lapply(seq_len(nrow(map)),function(i){
 values<-sapply(DONORS,function(d){a<-anns[[d]];paste(split_ids(a$KEGG_IDs[match(map[[paste0(d,'_feature_id')]][i],a$feature_id)]),collapse='/')})
 ref<-paste(split_ids(map$KEGG_IDs[i]),collapse='/');nonempty<-values[nzchar(values)]
 data.table(shared_id=map$shared_id[i],annotation=map$annotation[i],P39_KEGG=values[1],P41_KEGG=values[2],P59_KEGG=values[3],reference_KEGG=ref,nonempty_mappings_agree=length(unique(c(nonempty,if(nzchar(ref))ref)))<=1)
}));stopifnot(all(ma$nonempty_mappings_agree));fwrite(ma,'results/tables/KEGG_mapping_cross_donor_audit.csv')
ora<-function(a,universe){
 query<-intersect(unique(unlist(lapply(a$KEGG_IDs,split_ids))),universe);query<-sort(query)
 N<-length(universe);n<-length(query);stopifnot(n>0)
 d<-rbindlist(lapply(names(sets),function(id){
  bg<-intersect(unique(sets[[id]]),universe);M<-length(bg);if(M<2)return(NULL)
  hits<-intersect(bg,query);k<-length(hits)
  hf<-a[vapply(a$KEGG_IDs,function(x)any(split_ids(x)%in%hits),logical(1)),]
  data.table(ID=id,Description=pn$PathwayName[match(id,pn$PathwayID)],Count=k,GeneRatio=paste0(k,'/',n),BgRatio=paste0(M,'/',N),pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),FoldEnrich=(k/n)/(M/N),KEGG_IDs=paste(hits,collapse='/'),Feature_IDs=paste(hf$feature_id,collapse=';'))
 }))
 d[,p.adjust:=NA_real_];d[Count>0,p.adjust:=p.adjust(pvalue,'BH')];d[,p.adjust_all_eligible:=p.adjust(pvalue,'BH')]
 list(all=d,query=query,N=N,n=n)
}
# Numerically reproduce the supplied P59 panel_d before analyzing new lists.
s<-fread('inputs/P59_panel_d_reference/inputs/all_features_GAM_significance.csv');a<-copy(p59ann)
a[,q:=s$association_BH_FDR_Osi[match(feature_id,s$feature_id)]];a<-a[!is.na(q)&q<FDR]
r<-ora(a,universe);old<-fread('inputs/P59_panel_d_reference/tables/KEGG_enrichment_all_eligible_pathways.csv')
old<-old[match(r$all$ID,ID)];checks<-rbindlist(lapply(c('Count','pvalue','FoldEnrich','p.adjust','p.adjust_all_eligible'),function(col)data.table(metric=col,max_absolute_difference=max(abs(r$all[[col]]-old[[col]]),na.rm=TRUE))))
stopifnot(nrow(old)==nrow(r$all),all(checks$max_absolute_difference<1e-10));fwrite(checks,'results/validation/reproduce_P59_panel_d.csv')
queries<-list(at_least_two_donors=map[selected_at_least_two==TRUE])
summary<-list();combined<-list()
for(nm in names(queries)){
 q<-queries[[nm]];setorder(q,section_rank,pattern,row_order)
 q[,in_pathway_background:=vapply(KEGG_IDs,function(x)length(intersect(split_ids(x),universe))>0,logical(1))]
 q[,KEGG_in_background:=vapply(KEGG_IDs,function(x)paste(intersect(split_ids(x),universe),collapse='/'),character(1))]
 fwrite(q,paste0('results/tables/',nm,'_metabolite_list.csv'))
 r<-ora(q,universe);d<-r$all
 setorder(d,pvalue,-FoldEnrich,Description)
 fwrite(d,paste0('results/tables/',nm,'_KEGG_all_eligible.csv'))
 fwrite(data.table(KEGGID=r$query),paste0('results/tables/',nm,'_query_unique_KEGG_compounds.csv'))
 comp<-rbindlist(lapply(r$query,function(id){h<-q[vapply(KEGG_IDs,function(x)id%in%split_ids(x),logical(1))];data.table(KEGGID=id,shared_ids=paste(h$shared_id,collapse=';'),annotations=paste(h$annotation,collapse=';'),n_features=nrow(h))}))
 fwrite(comp,paste0('results/tables/',nm,'_compound_feature_mapping.csv'))
 display<-d[Count>0&pvalue<.05];setorder(display,p.adjust,-FoldEnrich,Description);display[,top_to_bottom_order:=seq_len(.N)]
 fwrite(display,paste0('results/tables/',nm,'_plot_data.csv'))
 lab<-if(nm=='three_donors')'Concordant and significant in all 3 donors' else 'All 3 donors concordant; at least 2 significant'
 if(nrow(display)){
  display[,Description:=factor(Description,levels=rev(Description))]
  p<-ggplot(display,aes(-log10(p.adjust),Description))+geom_segment(aes(x=0,xend=-log10(p.adjust),yend=Description),colour='#CCD5DF',linewidth=.7)+geom_point(aes(size=Count,fill=FoldEnrich),shape=21,colour='black',stroke=.2,alpha=.9)+geom_vline(xintercept=-log10(.05),linetype=2,colour='#B2182B')+scale_size_area(max_size=9,breaks=sort(unique(display$Count)))+scale_fill_gradientn(colours=c('#6AADC8','#F6D57E','#D1493F','#A50026'))+theme_classic(base_size=11)+theme(axis.title.y=element_blank(),plot.title.position='plot',plot.caption.position='plot',plot.title=element_text(face='bold'))+labs(x=expression(-log[10]('BH-FDR')),fill='Fold enrichment',size='Unique KEGG\ncompounds',title=paste('f |',lab),subtitle=sprintf('%d features -> %d unique KEGG compounds | full pathway background: %d compounds',nrow(q),r$n,r$N),caption='Shown: nominal P < 0.05. BH adjustment: hit pathways (P59 method). Dashed line: BH-FDR = 0.05.\nUpward and downward trends are pooled. Selection scope: 91 features shared by all 3 donors.')
 }else p<-ggplot()+annotate('text',x=0,y=0,label='No pathways with nominal P < 0.05')+theme_void()+labs(title=lab)
 h<-max(5.5,2.5+.34*nrow(display))
 ggsave(paste0('results/figures/panel_f_',nm,'_KEGG_enrichment.pdf'),p,width=13,height=h,device=cairo_pdf,bg='white')
 ggsave(paste0('results/figures/panel_f_',nm,'_KEGG_enrichment.png'),p,width=13,height=h,dpi=240,bg='white')
 ggsave(paste0('results/figures/panel_f_',nm,'_KEGG_enrichment.svg'),p,width=13,height=h,device=grDevices::svg,bg='white')
 summary[[nm]]<-data.table(list=nm,features=nrow(q),features_with_KEGG_in_background=sum(q$in_pathway_background),unique_query_compounds=r$n,background_compounds=r$N,eligible_pathways=nrow(d),hit_pathways=sum(d$Count>0),nominal_P_lt_005=sum(d$pvalue<.05),hit_pathway_BH_FDR_lt_005=sum(d$p.adjust<.05,na.rm=TRUE),all_eligible_BH_FDR_lt_005=sum(d$p.adjust_all_eligible<.05))
 combined[[nm]]<-cbind(data.table(list=nm),d)
}
fwrite(rbindlist(summary),'results/tables/enrichment_summary.csv');fwrite(rbindlist(combined),'results/tables/KEGG_results.csv')
fwrite(map[,.(shared_id,annotation,selected_three,selected_at_least_two,n_significant_donors,concordant,n_qualifying_pairs,qualifying_pairs,qualifying_shapes)],'results/tables/list_membership_all_91_features.csv')
writeLines(capture.output(sessionInfo()),'results/validation/R_sessionInfo.txt')
print(rbindlist(summary));for(nm in names(combined)){cat('\n',nm,'\n');print(combined[[nm]][p.adjust<.05,.(Description,Count,FoldEnrich,pvalue,p.adjust,p.adjust_all_eligible)])}
