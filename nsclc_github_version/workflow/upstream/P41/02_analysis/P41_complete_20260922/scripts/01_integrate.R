suppressPackageStartupMessages({library(data.table);library(uwot);library(openxlsx)})
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));d<-readRDS(file.path(root,'audit/raw_import.rds'))
out<-function(x,n)fwrite(x,file.path(root,'audit',n),na='')
for(p in c('a','b','c','d','e','g','h'))for(sub in c('inputs','scripts','tables','figures'))dir.create(file.path(root,paste0('panel_',p),sub),recursive=TRUE,showWarnings=FALSE)
# Validate actual scan IDs against updated source annotation.
annpath<-Sys.getenv('NSCLC_P41_ANNOTATION_CSV', unset=file.path(root,'inputs','original_annotation.csv'))
ann<-fread(annpath);m1<-d$metadata$B1;ix<-match(m1$scan_id,ann$scan_Id);expected<-sub('^Day','D',ann$Type[ix]);ok<-(is.na(m1$day)&is.na(expected))|(!is.na(m1$day)&!is.na(expected)&m1$day==expected)
out(data.table(cell_id=m1$cell_id,scan_id=m1$scan_id,h5ad_day=m1$day,source_day=expected,matched=!is.na(ix),label_equal=ok),'updated_annotation_alignment.csv');stopifnot(all(ok))
records<-copy(d$records);PPM_TOL<-5
setorder(records,mz,bi,source_column);n<-nrow(records);parent<-seq_len(n);members<-as.list(seq_len(n))
findroot<-function(i){while(parent[i]!=i){parent[i]<<-parent[parent[i]];i<-parent[i]};i}
edges<-vector('list',n)
for(i in seq_len(n-1L)){j<-i+1L;rr<-integer();while(j<=n&&records$mz[j]<=records$mz[i]*(1+PPM_TOL/1e6)){if(records$bi[j]!=records$bi[i])rr<-c(rr,j);j<-j+1L};if(length(rr))edges[[i]]<-data.table(ppm=2*(records$mz[rr]-records$mz[i])/(records$mz[rr]+records$mz[i])*1e6,i=i,j=rr)}
edges<-rbindlist(edges);setorder(edges,ppm,i,j);rejected<-list()
for(k in seq_len(nrow(edges))){i<-edges$i[k];j<-edges$j[k];u<-findroot(i);v<-findroot(j);if(u==v)next;rr<-c(members[[u]],members[[v]]);m<-records$mz[rr];if(anyDuplicated(records$bi[rr])||2*diff(range(m))/sum(range(m))*1e6>PPM_TOL){rejected[[length(rejected)+1L]]<-edges[k];next};parent[v]<-u;members[[u]]<-rr;members[v]<-list(NULL)}
groups<-members[lengths(members)>0];groups<-groups[order(vapply(groups,function(g)mean(records$mz[g]),numeric(1)))];out<-rbindlist(lapply(seq_along(groups),function(k){d<-copy(records[groups[[k]]]);d[,feature_id:=sprintf('F%05d',k)];d}))
mapping<-out
out<-function(x,n)fwrite(x,file.path(root,'audit',n),na='')
out(mapping,'feature_mapping.csv')
meta<-rbindlist(d$metadata)[excluded==FALSE];rawlist<-lapply(d$raw,function(x)x[rownames(x)%in%meta$cell_id,,drop=FALSE])
catalog<-mapping[,.(mz=mean(mz),n_series=.N,annotation_candidates=paste(sort(unique(na.omit(source_metabolite[nzchar(source_metabolite)]))),collapse='; ')),by=feature_id]
# Structural absence stays NA; observed nondetection is zero, supported by original n_metabolites.
rawfull<-matrix(NA_real_,nrow(meta),nrow(catalog),dimnames=list(meta$cell_id,catalog$feature_id))
for(b in names(rawlist)){mp<-mapping[batch==b];ii<-which(meta$batch==b);rawfull[ii,mp$feature_id]<-rawlist[[b]][meta$cell_id[ii],mp$source_column+1,drop=FALSE]}
rates<-sapply(unique(meta$group),function(g)colMeans(!is.na(rawfull[meta$group==g,,drop=FALSE])&rawfull[meta$group==g,,drop=FALSE]>0))
support<-sapply(unique(meta$batch),function(b)colSums(rawfull[meta$batch==b,,drop=FALSE]>0,na.rm=TRUE)>0)
catalog[,`:=`(n_support_series=rowSums(support),max_label_detection=apply(rates,1,max),all_three=n_series==3)]
catalog[,keep:=all_three|(n_support_series>=2&max_label_detection>=.5)];catalog[,actual_model_feature:=all_three]
out(catalog,'feature_catalog_screening.csv');out(data.table(feature_id=catalog$feature_id,rates),'feature_detection_by_label.csv')
# Retain supplementary catalog + NA-bearing raw matrix, never zero-fill unmeasured series.
fwrite(data.table(cell_id=meta$cell_id,rawfull[,catalog$keep,drop=FALSE]),file.path(root,'matrix/screened_catalog_raw_with_structural_NA.csv'),na='NA')
fids<-catalog[all_three==TRUE,feature_id];raw<-rawfull[,fids,drop=FALSE];stopifnot(!anyNA(raw))
X<-log1p(raw/meta$TIC*1e7);COMBAT_TOL<-1e-4
combat_affine<-function(x,batch,mask){
 lev<-sort(unique(batch[mask]));valid<-Reduce('&',lapply(lev,function(b)apply(x[mask&batch==b,,drop=FALSE],2,var)>1e-10))
 z<-x[mask,valid,drop=FALSE];bb<-batch[mask];center<-colMeans(z);fit<-z
 for(b in lev)fit[bb==b,]<-matrix(colMeans(z[bb==b,,drop=FALSE]),sum(bb==b),ncol(z),byrow=TRUE)
 vp<-colMeans((z-fit)^2);stopifnot(all(vp>0));s<-sweep(sweep(z,2,center,'-'),2,sqrt(vp),'/');out<-list()
 for(b in lev){
   sb<-s[bb==b,,drop=FALSE];n<-nrow(sb);gh<-colMeans(sb);dh<-apply(sb,2,var);dh[dh==0]<-1
   gb<-mean(gh);t2<-var(gh);dm<-mean(dh);dv<-var(dh);ap<-(2*dv+dm^2)/dv;bp<-(dm*dv+dm^3)/dv
   go<-gh;do<-dh;change<-Inf;iter<-0L
   while(change>COMBAT_TOL&&iter<1000L){gn<-(t2*n*gh+do*gb)/(t2*n+do);ss<-colSums(sweep(sb,2,gn,'-')^2);dn<-(.5*ss+bp)/(n/2+ap-1);change<-max(abs(gn-go)/go,abs(dn-do)/do,na.rm=TRUE);go<-gn;do<-dn;iter<-iter+1L}
   slope<-rep(1,ncol(x));offset<-rep(0,ncol(x));slope[valid]<-1/sqrt(dn);offset[valid]<-center*(1-slope[valid])-gn*sqrt(vp)*slope[valid]
   if(any(!valid)){targetmean<-Reduce('+',lapply(lev,function(l)colMeans(x[mask&batch==l,!valid,drop=FALSE])))/length(lev);offset[!valid]<-targetmean-colMeans(x[mask&batch==b,!valid,drop=FALSE])}
   out[[b]]<-list(slope=slope,offset=offset,iterations=iter)
 };out
}
apply_affine<-function(x,batch,transforms){out<-x;for(b in names(transforms)){ix<-batch==b;t<-transforms[[b]];out[ix,]<-sweep(sweep(x[ix,,drop=FALSE],2,t$slope,'*'),2,t$offset,'+')};out}
tr0<-combat_affine(X,meta$batch,meta$day=='D0');base<-apply_affine(X,meta$batch,tr0)
tr14<-combat_affine(base,meta$batch,meta$batch!='B3'&meta$day=='D14');second<-apply_affine(base,meta$batch,tr14)
w<-ifelse(meta$batch=='B3',0,pmin(meta$day_numeric/14,1));after<-base+(second-base)*w
stopifnot(all(is.finite(after)))
params<-rbindlist(lapply(c('D0','D14'),function(st){tt<-if(st=='D0')tr0 else tr14;rbindlist(lapply(names(tt),function(b)data.table(stage=st,batch=b,feature_id=fids,slope=tt[[b]]$slope,offset=tt[[b]]$offset)))}));out(params,'anchor_affine_parameters.csv')
set.seed(42);pc<-prcomp(X,center=TRUE,scale.=FALSE,rank.=50);emb<-function(z){set.seed(42);uwot::umap(z,n_neighbors=15,min_dist=.3,n_threads=1,n_sgd_threads=1,verbose=FALSE)}
u0<-emb(pc$x);u1<-emb(scale(after,center=pc$center,scale=FALSE)%*%pc$rotation);ud<-emb(scale(base,center=pc$center,scale=FALSE)%*%pc$rotation)
meta[,`:=`(UMAP_before_1=u0[,1],UMAP_before_2=u0[,2],UMAP_after_1=u1[,1],UMAP_after_2=u1[,2],UMAP_D0only_1=ud[,1],UMAP_D0only_2=ud[,2])]
res<-Sys.getenv('NSCLC_RESOURCE_DIR', unset=file.path(root,'inputs','annotation_resources'));rr<-list.files(res,pattern='^[^.].*(xlsx|txt)$',full.names=TRUE);file.copy(rr,file.path(root,'resources'),overwrite=FALSE)
met<-read.xlsx(file.path(root,'resources/Detected features and metabolites.xlsx'),sheet='Metabolites')
write.csv(met,file.path(root,'resources/metabolites.csv'),row.names=FALSE)
norm<-function(x)gsub('[^a-z0-9]','',tolower(iconv(x,to='ASCII//TRANSLIT',sub='')))
split_ids<-function(x){z<-unlist(strsplit(paste(na.omit(x),collapse='/'),'[/;]'));sort(unique(trimws(z[grepl('^C[0-9]{5}$',trimws(z))])))}
ft<-catalog[match(fids,feature_id)];ft[,`:=`(annotation=annotation_candidates,annotation_name=annotation_candidates,mz_center=mz)]
# Split source candidate strings, preserving all candidate identities and no name-based peak merging.
ft[,KEGG_IDs:=vapply(annotation,function(x){tokens<-trimws(unlist(strsplit(x,';| \\| ')));paste(split_ids(met$KEGG[norm(met$Metabolite.name)%in%norm(tokens)]),collapse='/')},'')]
pn<-read.xlsx(file.path(root,'resources/has_MetabolicPathway.xlsx'),sheet='hsa',colNames=FALSE)[,1:2];names(pn)<-c('PathwayID','PathwayName');pn$PathwayID<-sub('^path:','',pn$PathwayID)
path<-fread(file.path(root,'resources/path.cpd.txt'),header=FALSE);setnames(path,c('PathwayID','KEGGID'));path[,PathwayID:=sub('^path:','',PathwayID)];path[,KEGGID:=sub('^cpd:','',KEGGID)];path<-unique(path[PathwayID%in%pn$PathwayID&grepl('^C[0-9]{5}$',KEGGID)])
ft[,PathwayIDs:=vapply(KEGG_IDs,function(x)paste(sort(unique(path[KEGGID%in%split_ids(x),PathwayID])),collapse=';'),'')]
ft[,PathwayNames:=vapply(PathwayIDs,function(x)paste(na.omit(pn$PathwayName[match(strsplit(x,';',fixed=TRUE)[[1]],pn$PathwayID)]),collapse='; '),'')]
ft[,mapping_note:='Candidate name mapping to supplied Resource; generic lipid IDs do not establish structure; distinct m/z retained.']
fwrite(path,file.path(root,'panel_d/inputs/TERM2GENE.csv'));fwrite(pn,file.path(root,'panel_d/inputs/TERM2NAME.csv'))
fwrite(ft,file.path(root,'matrix/P41_feature_annotations.csv'));fwrite(meta,file.path(root,'panel_a/inputs/cells_and_UMAP_coordinates.csv'))
saveRDS(list(raw=raw,before=X,after=after,D0_only=base,metadata=meta,features=ft,feature_ids=fids,anchor_D0=tr0,anchor_D14=tr14),file.path(root,'audit/integrated.rds'),compress=FALSE);saveRDS(pc,file.path(root,'audit/PCA_before.rds'),compress=FALSE)
out(meta[,.N,by=.(batch,day,group)],'QC_retained_counts.csv')
print(c(cells=nrow(after),catalog=nrow(catalog),screened=sum(catalog$keep),model=ncol(after),annotated=sum(nzchar(ft$annotation)),mapped=sum(nzchar(ft$KEGG_IDs))))
