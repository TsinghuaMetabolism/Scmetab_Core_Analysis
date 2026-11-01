# Refit the SAME model as panel h, retaining its original full-family FDR values.
# All inputs are local CSV. This exports dense curves; it does not redefine tests.
N_KNOTS<-6;N_GRID<-160
args<-grep('^--file=',commandArgs(FALSE),value=TRUE)
base<-dirname(dirname(normalizePath(sub('^--file=','',args[1]))))
suppressPackageStartupMessages(library(mgcv))
m<-read.csv(file.path(base,'inputs/reference_pseudotime_treated_cells.csv'))
s<-read.csv(file.path(base,'inputs/all_features_GAM_significance.csv'))
x<-read.csv(file.path(base,'inputs/corrected_intensity_cells_by_features.csv'),check.names=FALSE)
stopifnot(identical(m$cell_id,x$cell_id),all(m$day!='D0'),identical(names(x)[-1],s$feature_id))
m$condition<-factor(m$condition,levels=c('Osi','TCD+Osi'))
lo<-max(tapply(m$pseudotime,m$condition,quantile,.05));hi<-min(tapply(m$pseudotime,m$condition,quantile,.95))
g<-seq(lo,hi,length.out=N_GRID)
nd<-expand.grid(pseudotime=g,condition=factor(c('Osi','TCD+Osi'),levels=levels(m$condition)))
one<-function(i){
  dat<-m;dat$y<-x[[i+1]]
  fit<-gam(y~condition+s(pseudotime,by=condition,k=N_KNOTS,bs='cr'),data=dat,family=gaussian(),method='REML')
  st<-summary(fit)$s.table
  p<-st[c('s(pseudotime):conditionOsi','s(pseudotime):conditionTCD+Osi'),'p-value']
  list(pred=as.numeric(predict(fit,nd)),p=p)
}
res<-parallel::mclapply(seq_len(nrow(s)),one,mc.cores=4)
stopifnot(!any(vapply(res,inherits,logical(1),'try-error')))
P<-do.call(rbind,lapply(res,`[[`,'p'))
check<-data.frame(feature_id=s$feature_id,p_Osi=P[,1],p_TCDOsi=P[,2],
  max_p_difference=pmax(abs(P[,1]-s$association_p_Osi),abs(P[,2]-s$association_p_TCDOsi)))
stopifnot(max(check$max_p_difference)<1e-7)
Y<-do.call(rbind,lapply(res,`[[`,'pred'))
for(j in 1:2){
  mat<-Y[,seq_len(N_GRID)+(j-1)*N_GRID,drop=FALSE];colnames(mat)<-sprintf('T%03d',seq_len(N_GRID))
  write.csv(data.frame(feature_id=s$feature_id,mat),file.path(base,paste0('tables/GAM_fitted_',c('Osi','TCDOsi')[j],'.csv')),row.names=FALSE)
}
write.csv(check,file.path(base,'tables/GAM_refit_validation.csv'),row.names=FALSE)
write.csv(data.frame(column=sprintf('T%03d',seq_len(N_GRID)),pseudotime=g),file.path(base,'tables/GAM_pseudotime_grid.csv'),row.names=FALSE)
cat('Refitted',nrow(s),'features. Maximum P difference:',max(check$max_p_difference),'\n')
