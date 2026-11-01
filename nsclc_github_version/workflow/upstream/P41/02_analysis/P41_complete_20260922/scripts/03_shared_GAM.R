suppressPackageStartupMessages({library(data.table);library(mgcv);library(MASS)})
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));d<-readRDS(file.path(root,'audit/analyzed.rds'))
# CSV is the reproducible model interface.
ix<-d$metadata$day!='D0';m<-as.data.frame(d$metadata[ix]);m$condition<-factor(ifelse(m$batch=='B3','TCD+Osi','Osi'),levels=c('Osi','TCD+Osi'));X<-d$after[ix,,drop=FALSE];f<-as.data.frame(d$features)
wrmat<-function(X,path)fwrite(data.table(cell_id=rownames(X),X),path)
for(p in c('c','h'))fwrite(m,file.path(root,paste0('panel_',p), 'inputs/reference_pseudotime_treated_cells.csv'))
wrmat(X,file.path(root,'panel_c/inputs/corrected_intensity_cells_by_features.csv'));fwrite(f,file.path(root,'panel_c/inputs/metabolite_feature_annotations.csv'))
# Re-read CSV, so p-values and dense predictions match independent redraw/refits.
m<-read.csv(file.path(root,'panel_c/inputs/reference_pseudotime_treated_cells.csv'));m$condition<-factor(m$condition,levels=c('Osi','TCD+Osi'));xx<-read.csv(file.path(root,'panel_c/inputs/corrected_intensity_cells_by_features.csv'),check.names=FALSE);X<-as.matrix(xx[,-1]);rownames(X)<-xx$cell_id
lo<-max(tapply(m$pseudotime,m$condition,quantile,.05));hi<-min(tapply(m$pseudotime,m$condition,quantile,.95));stopifnot(is.finite(lo),is.finite(hi),hi>lo)
grid<-seq(lo,hi,length.out=160);nd<-expand.grid(pseudotime=grid,condition=factor(c('Osi','TCD+Osi'),levels=levels(m$condition)))
wald<-function(C,b,V){est<-drop(C%*%b);S<-C%*%V%*%t(C);r<-qr(S)$rank;if(r==0)return(c(stat=NA,df=NA,p=NA));w<-drop(t(est)%*%ginv(S)%*%est);c(stat=w,df=r,p=pchisq(w,r,lower.tail=FALSE))}
one<-function(i){
 dat<-m;dat$y<-X[,i];fit<-gam(y~condition+s(pseudotime,by=condition,k=6,bs='cr'),data=dat,family=gaussian(),method='REML');sm<-summary(fit)$s.table
 g8<-seq(lo,hi,length.out=8);Lm<-predict(fit,data.frame(condition=factor('Osi',levels=levels(m$condition)),pseudotime=g8),type='lpmatrix');Lc<-predict(fit,data.frame(condition=factor('TCD+Osi',levels=levels(m$condition)),pseudotime=g8),type='lpmatrix');b<-coef(fit);V<-vcov(fit,unconditional=TRUE);pat<-wald(Lc-Lm,b,V);sem<-wald(Lm[8,,drop=FALSE]-Lm[1,,drop=FALSE],b,V);sec<-wald(Lc[8,,drop=FALSE]-Lc[1,,drop=FALSE],b,V)
 pr<-predict(fit,nd,se.fit=TRUE,unconditional=TRUE)
 st<-data.table(feature_id=colnames(X)[i],fit_ok=TRUE,association_p_Osi=sm['s(pseudotime):conditionOsi','p-value'],association_p_TCDOsi=sm['s(pseudotime):conditionTCD+Osi','p-value'],pattern_stat=pat['stat'],pattern_df=pat['df'],pattern_pvalue=pat['p'],startEnd_delta_Osi=drop((Lm[8,]-Lm[1,])%*%b),startEnd_p_Osi=sem['p'],startEnd_delta_TCDOsi=drop((Lc[8,]-Lc[1,])%*%b),startEnd_p_TCDOsi=sec['p'],converged=fit$converged,edf_Osi=sm[1,'edf'],edf_TCDOsi=sm[2,'edf'])
 list(stats=st,curve=as.numeric(pr$fit),se=as.numeric(pr$se.fit))
}
cat('Fitting',ncol(X),'features on',nrow(X),'treated cells; support',lo,hi,'\n');rr<-parallel::mclapply(seq_len(ncol(X)),one,mc.cores=4);stopifnot(!any(vapply(rr,inherits,logical(1),'try-error')))
s<-rbindlist(lapply(rr,`[[`,'stats'));pcols<-c('association_p_Osi','association_p_TCDOsi','pattern_pvalue','startEnd_p_Osi','startEnd_p_TCDOsi');qcols<-c('association_BH_FDR_Osi','association_BH_FDR_TCDOsi','pattern_BH_FDR','startEnd_BH_FDR_Osi','startEnd_BH_FDR_TCDOsi')
for(i in seq_along(pcols))s[,(qcols[i]):=p.adjust(get(pcols[i]),'BH',n=ncol(X))];s[,BH_family_size:=ncol(X)]
for(p in c('c','d','h'))fwrite(s,file.path(root,paste0('panel_',p),'inputs/all_features_GAM_significance.csv'))
Y<-do.call(rbind,lapply(rr,`[[`,'curve'));for(j in 1:2){mat<-Y[,seq_len(160)+(j-1)*160,drop=FALSE];colnames(mat)<-sprintf('T%03d',1:160);fwrite(data.table(feature_id=s$feature_id,mat),file.path(root,paste0('panel_c/tables/GAM_fitted_',c('Osi','TCDOsi')[j],'.csv')))}
fwrite(data.table(column=sprintf('T%03d',1:160),pseudotime=grid),file.path(root,'panel_c/tables/GAM_pseudotime_grid.csv'))
fwrite(data.table(parameter=c('n_features','n_Osi','n_TCDOsi','common_lo','common_hi','n_grid','D0_included'),value=c(ncol(X),sum(m$condition=='Osi'),sum(m$condition!='Osi'),lo,hi,160,0)),file.path(root,'panel_c/tables/model_parameters.csv'))
fa<-f[!is.na(f$annotation)&nzchar(f$annotation),];for(p in c('d','h'))fwrite(fa,file.path(root,paste0('panel_',p),'inputs/annotated_feature_KEGG_pathway_mapping.csv'))
ch<-f$feature_id[grepl('(^|;)map05231(;|$)',f$PathwayIDs)];wrmat(X[,ch,drop=FALSE],file.path(root,'panel_h/inputs/choline_pathway_corrected_intensity_cells_by_features.csv'))
fwrite(merge(f[f$feature_id%in%ch,],s,by='feature_id'),file.path(root,'panel_h/tables/P41_choline_pathway_Osi_vs_OsiTCD_significance.csv'))
ix<-d$metadata$batch!='B3';wrmat(d$after[ix,,drop=FALSE],file.path(root,'panel_e/inputs/corrected_intensity_cells_by_features.csv'));wrmat(1*(d$raw[ix,,drop=FALSE]>0),file.path(root,'panel_e/inputs/raw_detection_cells_by_features.csv'));fwrite(d$metadata[ix],file.path(root,'panel_e/inputs/Osi_cell_metadata.csv'));fwrite(f,file.path(root,'panel_e/inputs/metabolite_feature_annotations.csv'))
print(s[,lapply(.SD,function(x)sum(x<.05)),.SDcols=qcols])
