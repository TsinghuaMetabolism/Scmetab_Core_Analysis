suppressPackageStartupMessages({library(rhdf5);library(data.table)})
args<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',args[1]))))
decode<-function(v){if(is.list(v)){if(all(c('categories','codes')%in%names(v))){ca<-decode(v$categories);out<-rep(NA_character_,length(v$codes));ok<-v$codes>=0;out[ok]<-ca[v$codes[ok]+1];return(out)};if('values'%in%names(v)){out<-decode(v$values);if('mask'%in%names(v))out[as.logical(v$mask)]<-NA;return(out)};stop('Unsupported encoding')};as.vector(v)}
hc<-function(f,p)decode(h5read(f,p))
# Original B2/B3 H5AD are missing from this handoff; the 3-file assertion prevents partial import.
files<-list.files(file.path(root,'../../01_raw'),pattern='^[^.].*h5ad$',full.names=TRUE)
stopifnot(length(files)==3,grepl('20260909-Annotated',files[1]))
manifest<-data.table(batch=paste0('B',1:3),source=files,bytes=file.info(files)$size,modified=as.character(file.info(files)$mtime),md5=unname(tools::md5sum(files)))
raw<-list();meta<-list();features<-list()
for(i in 1:3){f<-files[i];b<-manifest$batch[i];x<-t(h5read(f,'X'));day<-sub('^Day','D',hc(f,'obs/cell_type'));id<-hc(f,'obs/CellNumber');tic<-as.numeric(hc(f,'obs/TIC'));tc<-as.numeric(hc(f,'obs/TIC_corrected'))
 if(anyNA(x)){stopifnot(all(rowSums(!is.na(x))==hc(f,"obs/n_metabolites")));x[is.na(x)]<-0}
 m<-data.table(cell_id=paste('P41',b,id,sep='__'),original_cell_id=id,sample_id=paste('P41',b,day,sep='_'),donor='P41',batch=b,batch_short=b,day=day,day_numeric=as.numeric(sub('D','',day)),scan_id=hc(f,'obs/scan_Id'),scan_start_time=hc(f,'obs/scan_start_time'),TIC=tic,TIC_corrected=tc)
 m[,treatment_series:=ifelse(batch=='B3','Osi+TCD','Osi')];m[,group:=ifelse(day=='D0','untreated_D0',paste0(ifelse(batch=='B3','OsiTCD_','Osi_'),day))];m[,treatment_status:=ifelse(day=='D0','untreated',treatment_series)]
 stopifnot(nrow(x)==nrow(m),!anyDuplicated(m$cell_id),all(is.finite(x)),all(x>=0),all(is.na(day)|day%in%c('D0','D3','D8','D14')))
 m[,`:=`(detected_features=rowSums(x>0),feature_sum=rowSums(x),top1_fraction=apply(x,1,max)/pmax(rowSums(x),1),log10_TIC=log10(TIC))]
 m[,excluded:=is.na(day)|!is.finite(TIC)|TIC<=0|!is.finite(TIC_corrected)|TIC_corrected<=0|detected_features==0]
 m[,exclusion_reason:=ifelse(excluded,'missing_day_or_nonpositive_or_invalid_TIC_corrected_or_empty_signal','retained')]
 for(v in c('log10_TIC','detected_features','top1_fraction'))m[,paste0('review_',v):={a<-get(v);med<-median(a,na.rm=TRUE);mad<-median(abs(a-med),na.rm=TRUE);if(mad==0)rep(FALSE,.N)else abs(a-med)>5*mad},by=day]
 rownames(x)<-m$cell_id;raw[[b]]<-x;meta[[b]]<-m
 features[[b]]<-data.table(batch=b,bi=i,source_column=seq_len(ncol(x))-1L,source_feature=hc(f,paste0('var/',h5readAttributes(f,'var')$`_index`)),mz=as.numeric(hc(f,'var/mz_center')),source_metabolite=hc(f,'var/metabolite'))
 manifest[i,`:=`(n_cells=nrow(x),n_features=ncol(x))]
}
fwrite(manifest,file.path(root,'audit/raw_source_manifest.csv'));mm<-rbindlist(meta);fwrite(mm,file.path(root,'audit/all_cell_QC_decisions.csv'));fwrite(mm[,.N,by=.(batch,scan_id,day,group,excluded)],file.path(root,'audit/scan_label_counts.csv'));print(mm[,.N,by=.(batch,day,excluded)])
saveRDS(list(raw=raw,metadata=meta,records=rbindlist(features),manifest=manifest),file.path(root,'audit/raw_import.rds'),compress=FALSE)
h5closeAll()
