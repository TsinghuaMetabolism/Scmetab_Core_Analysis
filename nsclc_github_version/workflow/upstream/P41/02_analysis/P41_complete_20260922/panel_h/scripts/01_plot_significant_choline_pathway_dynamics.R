# Pathway membership AND pseudotime association, using the original all-feature BH families.
PATHWAY_ID<-'map05231';FDR_CUTOFF<-.05
SIGNIFICANCE_RULE<-'either' # either, Osi, or pattern
N_KNOTS<-6;N_GRID<-160;COMMON_RANGE_QUANTILES<-c(.05,.95)
FACET_COLUMNS<-5;LINE_WIDTH<-.75;CI_ALPHA<-.15
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
read_input<-function(name)read.csv(file.path('inputs',name),check.names=FALSE,stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(mgcv)})
f<-read_input('annotated_feature_KEGG_pathway_mapping.csv');s<-read_input('all_features_GAM_significance.csv')
m<-read_input('reference_pseudotime_treated_cells.csv');inp<-read_input('choline_pathway_corrected_intensity_cells_by_features.csv')
Y<-as.matrix(inp[,-1,drop=FALSE]);rownames(Y)<-inp$cell_id
stopifnot(identical(rownames(Y),m$cell_id),!anyDuplicated(f$feature_id),all(m$day!='D0'),nrow(s)==length(unique(s$feature_id)))
f$in_choline_pathway<-vapply(strsplit(ifelse(is.na(f$PathwayIDs),'',f$PathwayIDs),';',fixed=TRUE),function(v)PATHWAY_ID%in%trimws(v),logical(1))
audit<-merge(f,s[,c('feature_id','association_BH_FDR_Osi','association_BH_FDR_TCDOsi','pattern_BH_FDR','association_p_Osi','association_p_TCDOsi')],by='feature_id',all.x=TRUE,sort=FALSE)
hit<-switch(SIGNIFICANCE_RULE,either=audit$association_BH_FDR_Osi<FDR_CUTOFF|audit$association_BH_FDR_TCDOsi<FDR_CUTOFF,Osi=audit$association_BH_FDR_Osi<FDR_CUTOFF,pattern=audit$pattern_BH_FDR<FDR_CUTOFF,stop('Unknown significance rule'))
audit$selected<-audit$in_choline_pathway & !is.na(hit) & hit
a<-audit[audit$selected,];a<-a[order(a$mz,a$feature_id),]
stopifnot(nrow(a)>0,all(a$feature_id%in%colnames(Y)))
# Preserve the supplied Resource mapping, including generic lipid IDs / ambiguous annotations.
a$mapping_note<-ifelse(a$KEGG_IDs=='C00157','Generic lipid KEGG mapping; preserve feature-level annotation uncertainty','Resource KEGG mapping')
write.csv(audit,'tables/pathway_significance_filter_audit.csv',row.names=FALSE,na='')
write.csv(a,'tables/selected_choline_pathway_features.csv',row.names=FALSE,na='')
m$condition<-factor(m$condition,levels=c('Osi','TCD+Osi'))
lo<-max(tapply(m$pseudotime,m$condition,quantile,COMMON_RANGE_QUANTILES[1]));hi<-min(tapply(m$pseudotime,m$condition,quantile,COMMON_RANGE_QUANTILES[2]))
nd<-expand.grid(pseudotime=seq(lo,hi,length.out=N_GRID),condition=factor(c('Osi','TCD+Osi'),levels=levels(m$condition)))
curves<-list();validation<-list()
for(i in seq_len(nrow(a))){
 id<-a$feature_id[i];dat<-m;dat$y<-Y[,id]
 fit<-gam(y~condition+s(pseudotime,by=condition,k=N_KNOTS,bs='cr'),data=dat,family=gaussian(),method='REML')
 pred<-predict(fit,nd,se.fit=TRUE,unconditional=TRUE)
 curves[[i]]<-data.frame(feature_id=id,annotation=a$annotation[i],mz=a$mz[i],KEGG_IDs=a$KEGG_IDs[i],nd,fit=as.numeric(pred$fit),lower=as.numeric(pred$fit-1.96*pred$se.fit),upper=as.numeric(pred$fit+1.96*pred$se.fit))
 st<-summary(fit)$s.table;pm<-st[grep('conditionOsi$',rownames(st)),'p-value'];pc<-st[grep('conditionTCD\\+Osi$',rownames(st)),'p-value']
 validation[[i]]<-data.frame(feature_id=id,max_association_p_difference=max(abs(c(pm-a$association_p_Osi[i],pc-a$association_p_TCDOsi[i]))))
}
cur<-do.call(rbind,curves);chk<-do.call(rbind,validation)
stopifnot(max(chk$max_association_p_difference)<1e-7)
write.csv(cur,'tables/pathway_GAM_curves_and_CI.csv',row.names=FALSE)
write.csv(chk,'tables/GAM_curve_refit_validation.csv',row.names=FALSE)
qtext<-function(q)ifelse(q<.001,'<0.001',paste0('=',formatC(q,format='f',digits=3)))
short<-ifelse(nchar(a$annotation)>53,paste0(substr(a$annotation,1,50),'...'),a$annotation)
wrap<-function(x){start<-seq(1,nchar(x),by=33);paste(substring(x,start,pmin(start+32,nchar(x))),collapse='\n')}
labels<-setNames(paste0(vapply(paste0(a$feature_id,' | ',short),wrap,character(1)),'\nOsi q',qtext(a$association_BH_FDR_Osi),' | TCD q',qtext(a$association_BH_FDR_TCDOsi)),a$feature_id)
cur$feature_id<-factor(cur$feature_id,levels=a$feature_id)
pal<-c(Osi='#2166AC','TCD+Osi'='#A23B72')
p<-ggplot(cur,aes(pseudotime,fit,colour=condition,fill=condition))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=CI_ALPHA,colour=NA)+geom_line(linewidth=LINE_WIDTH)+facet_wrap(~feature_id,ncol=FACET_COLUMNS,scales='free_y',labeller=labeller(feature_id=labels))+scale_colour_manual(values=pal)+scale_fill_manual(values=pal)+theme_classic(base_size=10)+theme(legend.position='top',strip.background=element_rect(fill='#EEF2F5',colour=NA),strip.text=element_text(size=8,lineheight=1.05),plot.title.position='plot',plot.caption.position='plot',plot.caption=element_text(hjust=0),plot.title=element_text(face='bold'),panel.spacing=grid::unit(1,'lines'))+labs(title='P41 | Significant dynamics within the choline pathway',subtitle=paste0('KEGG ',PATHWAY_ID,' | ',nrow(a),' annotated features | ',SIGNIFICANCE_RULE,' association significance rule; BH-FDR < ',FDR_CUTOFF),x='Osi-reference pseudotime',y='Corrected log-intensity (feature-specific scale)',colour='Treatment',fill='Treatment',caption='Fixed Osi reference; Osi+TCD uses reference-projected pseudotime. Gaussian GAM, k = 6; D0 excluded.\nFDR retains the original all-feature families. Bands: pointwise 95% intervals conditional on fixed pseudotime.\nMultiple lipid features map to generic KEGG C00157; feature IDs and full annotations are retained in CSV.')
ggsave('figures/panel_h_significant_choline_pathway_dynamics.pdf',p,width=20,height=max(7,ceiling(nrow(a)/FACET_COLUMNS)*3+1.5),device=cairo_pdf,bg='white')
ggsave('figures/panel_h_significant_choline_pathway_dynamics.png',p,width=20,height=max(7,ceiling(nrow(a)/FACET_COLUMNS)*3+1.5),dpi=300,bg='white')
cat('Pathway features selected:',nrow(a),'| maximum refit p difference:',max(chk$max_association_p_difference),'\n')
