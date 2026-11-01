# Resolve a single-State linear reference using experimental baseline ordering.
# Never read historical reference coordinates or p-values to choose orientation.
suppressPackageStartupMessages({library(monocle);library(Biobase);library(igraph);library(data.table)})
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
folder<-file.path(normalizePath(args[1],mustWork=TRUE),'refitted_trajectory')
stem<-'P59_pooled_top100_nc50_sig0.1_s20260918'
cds<-readRDS(file.path(folder,paste0(stem,'_cds.rds')))
native<-pData(cds)
linear<-length(unique(native$State))==1 && sum(degree(minSpanningTree(cds))>2)==0
if(!linear){
  fwrite(data.frame(applied=FALSE,reason='Not a single-State linear reference; retain declared root policy'),file.path(folder,'orientation_audit.csv'))
  quit(status=0)
}
day<-as.numeric(sub('^D','',native$day))
stopifnot(all(is.finite(day)),sum(day==0)>=5)
endpoint<-max(day);stopifnot(endpoint>0,sum(day==endpoint)>=5)
fits<-list(orderCells(cds,reverse=FALSE),orderCells(cds,reverse=TRUE))
audit<-rbindlist(lapply(seq_along(fits),function(i){
  pd<-pData(fits[[i]])[rownames(native),]
  data.frame(reverse=(i==2),baseline_day=0,endpoint_day=endpoint,
    baseline_median=median(pd$Pseudotime[day==0]),endpoint_median=median(pd$Pseudotime[day==endpoint]),
    baseline_n=sum(day==0),endpoint_n=sum(day==endpoint),pt_span=diff(range(pd$Pseudotime)))
}))
audit$baseline_endpoint_gap<-audit$endpoint_median-audit$baseline_median
audit$eligible<-audit$baseline_endpoint_gap>pmax(1e-10,.01*audit$pt_span)
stopifnot(sum(audit$eligible)==1)
chosen<-which(audit$eligible);fit<-fits[[chosen]];pd<-pData(fit)
audit$selected<-seq_len(nrow(audit))==chosen
audit$applied<-TRUE;audit$historical_reference_used_for_selection<-FALSE
fwrite(audit,file.path(folder,'orientation_audit.csv'))
tr<-fread(file.path(folder,paste0(stem,'_trajectory_cells.csv')),data.table=FALSE)
tr$native_refit_pseudotime<-tr$Pseudotime
tr$Pseudotime<-pd[tr$cell_id,'Pseudotime']
fwrite(tr,file.path(folder,'biologically_oriented_trajectory_cells.csv'))
saveRDS(fit,file.path(folder,'biologically_oriented_cds.rds'),compress=FALSE)
capture.output(sessionInfo(),file=file.path(folder,'R_sessionInfo.txt'))
print(audit)
