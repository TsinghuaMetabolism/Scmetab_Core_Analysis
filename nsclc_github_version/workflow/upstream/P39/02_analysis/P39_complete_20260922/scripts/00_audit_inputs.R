suppressPackageStartupMessages({library(data.table);library(rhdf5)})
args<-grep('^--file=',commandArgs(),value=T);root<-dirname(dirname(normalizePath(sub('^--file=','',args[1]))));setwd(root)
decode<-function(v){if(!is.list(v))return(as.vector(v));if('codes'%in%names(v)){cats<-decode(v$categories);z<-rep(NA_character_,length(v$codes));i<-v$codes>=0;z[i]<-cats[v$codes[i]+1];return(z)};if('values'%in%names(v)){z<-decode(v$values);if('mask'%in%names(v))z[as.logical(v$mask)]<-NA;return(z)};stop('encoding')}
hc<-function(f,p)decode(h5read(f,p))
files<-list.files(file.path(root,'../../01_raw'),pattern='h5ad$',full.names=T);raws<-list();metas<-list();vars<-list();aud<-list()
for(i in seq_along(files)){
 f<-files[i];b<-paste0('B',i);o<-h5read(f,'obs');v<-h5read(f,'var');cat(b,'obs:',names(o),'\nvar:',names(v),'\n');x<-t(h5read(f,'X'));m<-as.data.table(lapply(o,decode));vv<-as.data.table(lapply(v,decode));print(dim(x));print(table(m$cell_type));print(head(vv,2));
 m[,`:=`(batch=b,original_cell_id=as.character(CellNumber),cell_id=paste0('P39_',b,'__',CellNumber),day=sub('^Day','D',cell_type),donor='P39')];m[,`:=`(sample_id=paste('P39',batch,day,sep='_'),group=ifelse(day=='D0','untreated_D0',paste0('Osi_',day)),treatment_status=ifelse(day=='D0','untreated','Osimertinib'),detected_features=rowSums(x>0,na.rm=T),feature_sum=rowSums(x,na.rm=T),top1_fraction=apply(x,1,max,na.rm=T)/pmax(rowSums(x,na.rm=T),1))]
 stopifnot(all(rowSums(x>0,na.rm=TRUE)==m$n_metabolites)); m[,raw_NaN_count:=rowSums(is.nan(x))]; x[is.nan(x)]<-0

 m[,invalid_raw:=rowSums(!is.finite(x)|x<0)>0];m[,exclude_low_signal:=!is.finite(TIC)|TIC<=0|detected_features==0|(!is.na(TIC_corrected)&TIC_corrected<=0)];m[,excluded:=invalid_raw|exclude_low_signal];m[,exclusion_reason:=ifelse(invalid_raw,'invalid raw intensity',ifelse(exclude_low_signal,'nonpositive original/corrected TIC or zero detected features','retained; no sufficient exclusion evidence'))]
 m[,exclude_missing_cell_type:=is.na(cell_type)|!nzchar(trimws(cell_type))];m[exclude_missing_cell_type==TRUE,`:=`(excluded=TRUE,exclusion_reason='User requested exclusion: missing cell_type')]
 m[,log10_TIC:=log10(pmax(TIC,1))];for(n in c('log10_TIC','detected_features','top1_fraction'))m[,paste0('review_',n):={z<-get(n);s<-median(abs(z-median(z)));abs(z-median(z))>5*s},by=day]
 m[,boundary_review:=rank(scan_start_time,ties.method='first')<=ceiling(.05*.N)|rank(-scan_start_time,ties.method='first')<=ceiling(.05*.N),by=day]
 stopifnot(!anyDuplicated(m$cell_id));rownames(x)<-m$cell_id
 aud[[i]]<-data.table(batch=b,file=f,md5=unname(tools::md5sum(f)),cells=nrow(x),features=ncol(x),raw_missing=sum(m$raw_NaN_count),raw_negative=sum(x<0,na.rm=T),duplicate_cell_ids=anyDuplicated(m$CellNumber),duplicate_scan_ids=anyDuplicated(m$scan_Id),low_signal=sum(m$exclude_low_signal),excluded=sum(m$excluded))
 raws[[b]]<-x;metas[[b]]<-m;vars[[b]]<-vv
}
fwrite(rbindlist(aud),'tables/raw_source_manifest.csv');fwrite(rbindlist(metas,fill=TRUE),'tables/all_cell_QC_decisions.csv');fwrite(rbindlist(metas,fill=TRUE)[,.N,by=.(batch,day,excluded)],'tables/QC_counts.csv');saveRDS(list(raw=raws,meta=metas,var=vars), 'data/audited_raw.rds',compress=F)
print(rbindlist(aud));print(rbindlist(metas,fill=TRUE)[,.(n=.N,TIC_min=min(TIC),TICcorrected_min=min(TIC_corrected),detected_min=min(detected_features),scan_min=min(scan_start_time),scan_max=max(scan_start_time)),by=.(batch,day)])
