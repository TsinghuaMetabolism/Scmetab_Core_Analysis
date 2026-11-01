suppressPackageStartupMessages({library(data.table);library(uwot);library(cluster);library(openxlsx)})
a<-grep('^--file=',commandArgs(),value=T);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));setwd(root)
d<-readRDS('data/audited_raw.rds');m<-rbindlist(d$meta,fill=T)[excluded==FALSE];m[is.na(day),`:=`(day='Unknown',group='Unknown',sample_id='P39_B2_Unknown',treatment_status='unknown')]
r<-rbindlist(lapply(names(d$var),function(b){v<-d$var[[b]];data.table(batch=b,source_column=seq_len(nrow(v)),mz=v$mz_center,source_feature=v$`_index`,annotation=v$metabolite,original_annotation=v$metabolite_original,source_KEGG=v$metabolite_kegg)}));setorder(r,mz,batch,source_column)
# Minimum-ppm greedy one-to-one matching, no same-batch merging; identical P59 5ppm criterion for two batches.
ed<-rbindlist(lapply(which(r$batch=='B1'),function(i){j<-which(r$batch=='B2'&abs(r$mz-r$mz[i])/r$mz[i]*1e6<=5.001);if(!length(j))return(NULL);z<-data.table(i=i,j=j,ppm=2*abs(r$mz[j]-r$mz[i])/(r$mz[j]+r$mz[i])*1e6);z[ppm<=5]}));setorder(ed,ppm,i,j);used<-rep(F,nrow(r));groups<-list()
for(k in seq_len(nrow(ed)))if(!used[ed$i[k]]&&!used[ed$j[k]]){ii<-c(ed$i[k],ed$j[k]);groups[[length(groups)+1]]<-ii;used[ii]<-T};groups<-c(groups,as.list(which(!used)));groups<-groups[order(sapply(groups,function(ii)mean(r$mz[ii])))];mp<-rbindlist(lapply(seq_along(groups),function(k){z<-copy(r[groups[[k]]]);z[,feature_id:=sprintf('P39_F%05d',k)];z}));fwrite(mp,'tables/feature_mapping.csv')
f<-mp[,.(mz=mean(mz),n_series=.N,annotation=paste(unique(annotation[!is.na(annotation)&nzchar(annotation)]),collapse=' / '),original_annotations=paste(unique(original_annotation),collapse=' | '),source_KEGG=paste(unique(source_KEGG),collapse=' / ')),by=feature_id]
f[,keep:=n_series==2];f[,reason:=ifelse(keep,'shared by both batches','single-series feature cannot meet two-series support')];fwrite(f,'tables/feature_catalog.csv');f<-f[keep==TRUE];f[,annotation_name:=annotation]
raw<-matrix(NA_real_,nrow(m),nrow(f),dimnames=list(m$cell_id,f$feature_id));for(b in names(d$raw)){ii<-which(m$batch==b);mpp<-mp[batch==b][match(f$feature_id,feature_id)];raw[ii,]<-d$raw[[b]][m$cell_id[ii],mpp$source_column,drop=F]};stopifnot(all(is.finite(raw)))
X<-log1p(raw/m$TIC*1e7)
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
COMBAT_TOL<-1e-4
tr0<-combat_affine(X,m$batch,m$day=='D0');base<-apply_affine(X,m$batch,tr0);tr14<-combat_affine(base,m$batch,m$day=='D14');second<-apply_affine(base,m$batch,tr14);w<-pmin(suppressWarnings(as.numeric(sub('D','',m$day)))/14,1);w[is.na(w)]<-0;dual<-base+(second-base)*w
pc<-prcomp(X,center=T,scale.=F,rank.=50);proj<-function(z)sweep(z,2,pc$center,'-')%*%pc$rotation; mats<-list(Before=X,D0_only=base,D0_D14=dual);pcs<-lapply(mats,proj)
knn<-function(z,k=15){dd<-as.matrix(dist(z));diag(dd)<-Inf;t(apply(dd,1,function(x)order(x)[seq_len(min(k,nrow(z)-1))]))}
metrics<-list();set.seed(42)
for(day in c('D0','D14'))for(rep in 1:5){ix<-which(m$day==day);n<-min(table(m$batch[ix]));take<-unlist(lapply(c('B1','B2'),function(b)sample(ix[m$batch[ix]==b],n)));bb<-m$batch[take];for(nm in names(pcs)){z<-pcs[[nm]][take,,drop=F];nn<-knn(z);p<-rowMeans(matrix(bb[nn],nrow(nn))=='B1');metrics[[length(metrics)+1]]<-data.table(anchor=day,repeat_id=rep,method=nm,batch_mixing=mean(1/(p^2+(1-p)^2)),batch_ASW=mean(abs(silhouette(as.integer(factor(bb)),dist(z))[,3])),centroid_gap=sqrt(sum((colMeans(z[bb=='B1',])-colMeans(z[bb=='B2',]))^2)))}}
met<-rbindlist(metrics);fwrite(met,'tables/integration_subsampling.csv');print(met[,lapply(.SD,mean),by=.(anchor,method),.SDcols=c('batch_mixing','batch_ASW','centroid_gap')])
bio<-rbindlist(lapply(c('D0_only','D0_D14'),function(nm) rbindlist(lapply(c('B1','B2'),function(b){ii<-which(m$batch==b&m$day!='Unknown');z<-pcs$Before[ii,];v<-pcs[[nm]][ii,];nn<-knn(z);nn2<-knn(v);jac<-sapply(seq_along(ii),function(i)length(intersect(nn[i,],nn2[i,]))/length(union(nn[i,],nn2[i,])));data.table(method=nm,batch=b,distance_spearman=cor(as.vector(dist(z)),as.vector(dist(v)),method='spearman'),knn_jaccard=mean(jac))}))));fwrite(bio,'tables/biological_conservation.csv');print(bio)
saveRDS(list(raw=raw,before=X,D0_only=base,D0_D14=dual,meta=m,features=f,pc=pc,transforms=list(D0=tr0,D14=tr14)), 'data/integration_candidates.rds',compress=F)
fwrite(met[,lapply(.SD,mean),by=.(anchor,method),.SDcols=c('batch_mixing','batch_ASW','centroid_gap')],'tables/integration_metrics.csv')
