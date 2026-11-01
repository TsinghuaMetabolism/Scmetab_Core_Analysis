src<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);setwd(dirname(dirname(normalizePath(src))))
suppressPackageStartupMessages(library(mgcv))
a<-read.csv('inputs/feature_annotations.csv');m<-read.csv('inputs/trajectory_cells.csv')
x<-read.csv('inputs/P59_corrected_cell_feature_matrix.csv',check.names=FALSE)
meta<-read.csv('inputs/P59_cell_annotations.csv');v<-read.csv('inputs/cell_intensity_long.csv')
stopifnot(ncol(x)==1929,nrow(x)==1645,!anyDuplicated(x$cell_id),!anyDuplicated(m$cell_id),setequal(x$cell_id,m$cell_id))
j<-match(m$cell_id,x$cell_id);Y<-as.matrix(x[j,-1]);ids<-colnames(Y)
mi<-match(m$cell_id,meta$cell_id)
stopifnot(identical(m$sample_id,meta$sample_id[mi]),identical(m$day,meta$day[mi]),all(is.finite(Y)))
delta<-Y[cbind(match(v$cell_id,m$cell_id),match(v$feature_id,ids))]-v$intensity
stopifnot(max(abs(delta))<1e-10)
m$condition<-factor(ifelse(m$cohort=='Osi_reference','Osi','Osi+TCD'),levels=c('Osi','Osi+TCD'))
rr<-range(m$pseudotime[m$condition=='Osi']);m$t<-(m$pseudotime-rr[1])/diff(rr)
write.csv(m,'tables/model_cell_metadata.csv',row.names=FALSE)
write.csv(data.frame(cells=nrow(Y),features=ncol(Y),selected_features=nrow(a),maximum_intensity_difference=max(abs(delta))),
          'tables/input_validation.csv',row.names=FALSE)
curve_test<-function(f,C){
  b<-drop(C%*%coef(f));S<-C%*%vcov(f,unconditional=TRUE)%*%t(C)
  e<-eigen(S,symmetric=TRUE);keep<-e$values>max(e$values)*1e-8
  w<-sum(drop(crossprod(e$vectors[,keep,drop=FALSE],b))^2/e$values[keep])
  c(Wald=w,df=sum(keep),p=pchisq(w,sum(keep),lower.tail=FALSE))
}
allres<-list();allcur<-list();warnings_log<-list()
for(scenario in c('include_D0','exclude_D0')){
  keep<-if(scenario=='include_D0')rep(TRUE,nrow(m)) else m$day!='D0'
  d<-m[keep,];ym<-Y[keep,,drop=FALSE]
  lo<-max(tapply(d$t,d$condition,quantile,.05));hi<-min(tapply(d$t,d$condition,quantile,.95))
  # Use only the overlapping central support for the between-curve contrast.
  grid<-seq(lo,hi,length.out=40)
  nd0<-data.frame(t=grid,condition=factor('Osi',levels=levels(m$condition)))
  nd1<-data.frame(t=grid,condition=factor('Osi+TCD',levels=levels(m$condition)))
  res<-vector('list',length(ids))
  for(i in seq_along(ids)){
    d$y<-ym[,i];ww<-character()
    f<-withCallingHandlers(gam(y~condition+s(t,by=condition,k=6,bs='cr'),data=d,method='REML',family=gaussian()),
         warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')})
    if(length(ww))warnings_log[[length(warnings_log)+1]]<-data.frame(scenario,feature_id=ids[i],warning=paste(unique(ww),collapse=' | '))
    stopifnot(isTRUE(f$converged),all(is.finite(coef(f))))
    sm<-summary(f)$s.table
    pa<-sm[grep('conditionOsi$',rownames(sm)),'p-value'];pb<-sm[grep('conditionOsi\\+TCD$',rownames(sm)),'p-value']
    C<-predict(f,nd1,type='lpmatrix')-predict(f,nd0,type='lpmatrix')
    ov<-curve_test(f,C);sh<-curve_test(f,sweep(C,2,colMeans(C)))
    res[[i]]<-data.frame(feature_id=ids[i],scenario,n_Osi=sum(d$condition=='Osi'),n_combination=sum(d$condition=='Osi+TCD'),
      association_p_Osi=pa,association_p_TCDOsi=pb,overall_p=ov['p'],shape_p=sh['p'],
      overall_Wald=ov['Wald'],overall_df=ov['df'],shape_df=sh['df'],
      EDF_Osi=sm[1,'edf'],EDF_TCDOsi=sm[2,'edf'],comparison_t_min=lo,comparison_t_max=hi,converged=f$converged)
    if(ids[i]%in%a$feature_id){
      nd<-do.call(rbind,lapply(levels(d$condition),function(g){r<-range(d$t[d$condition==g]);data.frame(t=seq(r[1],r[2],length.out=250),condition=g)}))
      nd$condition<-factor(nd$condition,levels=levels(d$condition));p<-predict(f,nd,se.fit=TRUE,unconditional=TRUE)
      allcur[[length(allcur)+1]]<-data.frame(feature_id=ids[i],scenario,nd,fit=as.numeric(p$fit),se=as.numeric(p$se.fit),lower=as.numeric(p$fit-1.96*p$se.fit),upper=as.numeric(p$fit+1.96*p$se.fit))
      saveRDS(f,paste0('tables/',scenario,'_',ids[i],'_model.rds'))
    }
    if(i%%100==0)cat(scenario,i,'/',length(ids),'\n')
  }
  res<-do.call(rbind,res)
  for(pair in list(c('association_p_Osi','association_BH_FDR_Osi'),c('association_p_TCDOsi','association_BH_FDR_TCDOsi'),c('overall_p','overall_BH_FDR'),c('shape_p','shape_BH_FDR'))){
    stopifnot(all(is.finite(res[[pair[1]]])),all(res[[pair[1]]]>=0),all(res[[pair[1]]]<=1))
    # Preserve raw numerical zeros, but use an explicit conservative numerical floor before BH.
    res[[pair[2]]]<-p.adjust(pmax(res[[pair[1]]],.Machine$double.eps),'BH')
  }
  res$BH_family_size<-1928
  write.csv(res,paste0('tables/all_1928_GAM_',scenario,'.csv'),row.names=FALSE)
  allres[[scenario]]<-res
}
all<-do.call(rbind,allres);cur<-do.call(rbind,allcur)
write.csv(cur,'tables/fitted_curves_95CI_both_scenarios.csv',row.names=FALSE)
selected<-merge(a,all,by='feature_id',sort=FALSE);selected<-selected[order(match(selected$feature_id,a$feature_id),selected$scenario),]
selected$pdf_filename<-paste0(gsub('/',' or ',selected$annotation,fixed=TRUE),'.pdf')
write.csv(selected,'tables/selected_25_GAM_both_scenarios.csv',row.names=FALSE)
st<-selected[selected$scenario=='include_D0',];write.csv(st,'tables/curve_difference_statistics.csv',row.names=FALSE)
if(length(warnings_log))write.csv(do.call(rbind,warnings_log),'tables/model_warnings.csv',row.names=FALSE)
capture.output(sessionInfo(),file='tables/R_sessionInfo.txt')
cat('FINISHED\n');print(aggregate(cbind(association_BH_FDR_Osi,association_BH_FDR_TCDOsi,overall_BH_FDR)<.05~scenario,selected,sum))
