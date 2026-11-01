src<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(src))))
suppressPackageStartupMessages(library(mgcv))
sel<-read.csv('tables/selected_25_GAM_both_scenarios.csv');cur<-read.csv('tables/fitted_curves_95CI_both_scenarios.csv')
checks<-list()
for(sc in c('include_D0','exclude_D0')){
 all<-read.csv(paste0('tables/all_1928_GAM_',sc,'.csv'))
 stopifnot(nrow(all)==1928,!anyDuplicated(all$feature_id),all(all$converged),all(all$BH_family_size==1928))
 for(pair in list(c('association_p_Osi','association_BH_FDR_Osi'),c('association_p_TCDOsi','association_BH_FDR_TCDOsi'),c('overall_p','overall_BH_FDR'),c('shape_p','shape_BH_FDR')))
  stopifnot(max(abs(p.adjust(pmax(all[[pair[1]]],.Machine$double.eps),'BH')-all[[pair[2]]]))<1e-12)
 s<-sel[sel$scenario==sc,];stopifnot(nrow(s)==25)
 for(i in seq_len(nrow(s))){
  ix<-match(s$feature_id[i],all$feature_id)
  stopifnot(abs(s$association_BH_FDR_Osi[i]-all$association_BH_FDR_Osi[ix])<1e-12,
            abs(s$association_BH_FDR_TCDOsi[i]-all$association_BH_FDR_TCDOsi[ix])<1e-12)
  f<-readRDS(paste0('tables/',sc,'_',s$feature_id[i],'_model.rds'))
  cc<-cur[cur$scenario==sc & cur$feature_id==s$feature_id[i],]
  pr<-predict(f,cc,se.fit=TRUE,unconditional=TRUE)
  err<-max(abs(pr$fit-cc$fit),abs(pr$fit-1.96*pr$se.fit-cc$lower),abs(pr$fit+1.96*pr$se.fit-cc$upper))
  stopifnot(err<1e-10)
  sm<-summary(f)$s.table
  stopifnot(abs(sm[1,'p-value']-s$association_p_Osi[i])<1e-12,abs(sm[2,'p-value']-s$association_p_TCDOsi[i])<1e-12)
  checks[[length(checks)+1]]<-data.frame(scenario=sc,feature_id=s$feature_id[i],maximum_prediction_error=err,BH_full_family_verified=TRUE,association_p_verified=TRUE)
 }
}
write.csv(do.call(rbind,checks),'tables/validation.csv',row.names=FALSE)
cat('Verified both 1928-feature families, all 50 models, association P values, predictions and confidence intervals.\n')
