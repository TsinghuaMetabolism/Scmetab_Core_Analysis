# Translation to R of pseudobulk_de: cell-level one-day-vs-rest MWU.
FDR_CUTOFF<-.05;MIN_DETECTED_CELLS<-30L;MIN_SD<-1e-8;Z_LIMIT<-2
MIN_ABS_RANK_BISERIAL<-.30
DAYS<-c('D0','D3','D8','D14')
source_file<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(source_file))))
read_input<-function(name)read.csv(file.path('inputs',name),check.names=FALSE,stringsAsFactors=FALSE)
read_matrix<-function(name){x<-read_input(name);ids<-x[[1]];m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-ids;storage.mode(m)<-'double';m}
suppressPackageStartupMessages({library(pheatmap);library(grid);library(ragg);library(seriation)})
d<-list(cells=read_input('Osi_cell_metadata.csv'),corrected=read_matrix('corrected_intensity_cells_by_features.csv'),detected=read_matrix('raw_detection_cells_by_features.csv')>0,features=read_input('metabolite_feature_annotations.csv'));X<-d$corrected;m<-d$cells;f<-d$features;det<-d$detected
stopifnot(identical(rownames(X),m$cell_id),identical(dimnames(X),dimnames(det)),!anyDuplicated(m$cell_id))
f<-f[match(colnames(X),f$feature_id),];f$n_detected<-colSums(det);f$corrected_sd<-apply(X,2,sd)
# All input features are common to all batches: none structurally absent.
f$eligible<-f$n_detected>=MIN_DETECTED_CELLS&f$corrected_sd>MIN_SD
write.csv(f,'tables/feature_screen.csv',row.names=FALSE,na='')
X<-X[,f$eligible,drop=FALSE];det<-det[,f$eligible,drop=FALSE];f<-f[f$eligible,]
N<-nrow(X);ranks<-apply(X,2,rank,ties.method='average')
ties<-apply(X,2,function(x){t<-as.numeric(table(x));sum(t^3-t)})
tests<-do.call(rbind,lapply(DAYS,function(day){
 s<-m$day==day;n1<-sum(s);n2<-sum(!s);U<-colSums(ranks[s,,drop=FALSE])-n1*(n1+1)/2
 variance<-n1*n2/12*(N+1-ties/(N*(N-1)))
 z<-(abs(U-n1*n2/2)-.5)/sqrt(variance);p<-pmin(1,2*pnorm(z,lower.tail=FALSE));p[variance==0]<-1
 data.frame(feature_id=colnames(X),day=day,n_target=n1,n_rest=n2,U=U,p_value=p,rank_biserial=2*U/(n1*n2)-1,mean_target=colMeans(X[s,,drop=FALSE]),mean_rest=colMeans(X[!s,,drop=FALSE]),mean_difference=colMeans(X[s,,drop=FALSE])-colMeans(X[!s,,drop=FALSE]),detection_target=colMeans(det[s,,drop=FALSE]),detection_rest=colMeans(det[!s,,drop=FALSE]),row.names=NULL)
}))
tests$BH_FDR<-p.adjust(tests$p_value,'BH');tests$significant<-tests$BH_FDR<FDR_CUTOFF
tests$passes_FDR_and_effect<-tests$significant & abs(tests$rank_biserial)>=MIN_ABS_RANK_BISERIAL
# Independent check against R's tie-corrected, continuity-corrected asymptotic Wilcoxon.
check_ids<-seq_len(min(12,ncol(X)))
maxerr<-0
for(j in check_ids)for(day in DAYS){s<-m$day==day;pv<-wilcox.test(X[s,j],X[!s,j],exact=FALSE,correct=TRUE)$p.value;vv<-tests$p_value[tests$feature_id==colnames(X)[j]&tests$day==day];maxerr<-max(maxerr,abs(pv-vv))}
stopifnot(maxerr<1e-10)
write.csv(tests,'tables/all_contrasts.csv',row.names=FALSE)
M<-sapply(DAYS,function(day)colMeans(X[m$day==day,,drop=FALSE]));rownames(M)<-colnames(X)
f$min_contrast_BH_FDR<-vapply(f$feature_id,function(id)min(tests$BH_FDR[tests$feature_id==id]),numeric(1))
f$significant<-f$min_contrast_BH_FDR<FDR_CUTOFF
f$selected_FDR_and_effect<-f$feature_id %in% tests$feature_id[tests$passes_FDR_and_effect]
f$qualifying_days<-vapply(f$feature_id,function(id)paste(tests$day[tests$feature_id==id & tests$passes_FDR_and_effect],collapse=';'),character(1))
f$max_abs_rank_biserial<-vapply(f$feature_id,function(id)max(abs(tests$rank_biserial[tests$feature_id==id])),numeric(1))
for(day in DAYS)f[[paste0('mean_',day)]]<-M[,day]
write.csv(f,'tables/feature_summary.csv',row.names=FALSE,na='')
keep<-f$selected_FDR_and_effect&!is.na(f$annotation_name)&nzchar(trimws(f$annotation_name));a<-f[keep,];Y<-M[keep,,drop=FALSE]
stopifnot(nrow(Y)>=2)
sdpop<-sqrt(rowMeans((Y-rowMeans(Y))^2));Z<-(Y-rowMeans(Y))/sdpop
dst<-dist(Z);hc<-hclust(dst,method='ward.D2')
ord<-seriate(dst,method='OLO',control=list(hclust=hc));tree<-ord[[1]];stopifnot(inherits(tree,'hclust'))
idx<-get_order(ord);a<-a[idx,];zp<-Z[idx,,drop=FALSE]
write.csv(data.frame(a,row_order=seq_len(nrow(a)),setNames(as.data.frame(zp),paste0('rowZ_',DAYS)),check.names=FALSE),'tables/heatmap_data.csv',row.names=FALSE,na='')
write.csv(m,'tables/cells.csv',row.names=FALSE)
write.csv(data.frame(metric=c('cells','eligible_features','BH_family_size','FDR_significant_all_features','FDR_significant_annotated_features','FDR_and_effect_all_features','displayed_annotated_features','min_abs_rank_biserial','MWU_formula_vs_wilcox_max_difference'),value=c(nrow(X),ncol(X),nrow(tests),sum(f$significant),sum(f$significant&!is.na(f$annotation_name)&nzchar(f$annotation_name)),sum(f$selected_FDR_and_effect),nrow(a),MIN_ABS_RANK_BISERIAL,maxerr)),'tables/analysis_summary.csv',row.names=FALSE)
labels<-ifelse(nchar(f$annotation_name[keep])>75,paste0(substr(f$annotation_name[keep],1,72),'...'),f$annotation_name[keep])
draw<-function()pheatmap(pmax(pmin(Z,Z_LIMIT),-Z_LIMIT),cluster_rows=tree,cluster_cols=FALSE,labels_row=labels,labels_col=DAYS,border_color=NA,fontsize_row=9,fontsize_col=11,color=colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(100),breaks=seq(-Z_LIMIT,Z_LIMIT,length.out=101),main=paste0('P41 | Timepoint vs rest MWU\n',nrow(a),' annotated features | BH < ',FDR_CUTOFF,' and |rank-biserial| >= ',MIN_ABS_RANK_BISERIAL))
w<-13;h<-max(7,nrow(a)*.22+2)
cairo_pdf('figures/panel_e_pseudobulk_heatmap.pdf',width=w,height=h);draw();dev.off();agg_png('figures/panel_e_pseudobulk_heatmap.png',width=w,height=h,units='in',res=260);draw();dev.off()
print(read.csv('tables/analysis_summary.csv'))
