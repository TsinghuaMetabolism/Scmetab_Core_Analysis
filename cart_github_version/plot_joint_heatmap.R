# Cross-donor annotated dynamics on frozen per-donor curves.
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
source(file.path(dirname(normalizePath(script)),'plot_common.R'))
settings <- fromJSON(args[2])

dd <- file.path(panel_dir,'OriginalData');qa <- file.path(dd,'qa');dir.create(qa,showWarnings=FALSE)
r <- read.delim(file.path(dd,'row_metadata.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
d <- read.delim(file.path(dd,'heatmap_data.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
donors <- settings$reference_donors;titles <- settings$reference_labels
stopifnot(length(donors)==3,length(titles)==3)
groups <- unique(r$group);rowh <- 5.6;xs <- c(101,155,209);hw <- 44
width_mm <- 270;height_mm <- 14+length(groups)*10+nrow(r)*rowh+27
pal <- colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(257)
colors <- function(z)pal[1L+round((pmax(-2,pmin(2,z))+2)/4*256)]
wrap_name <- function(s,width=92,size=7) {
  s <- gsub('α','alpha',gsub('γ','gamma',s,fixed=TRUE),fixed=TRUE)
  gp <- gpar(fontsize=size,fontfamily='Helvetica');ans <- character();left <- s
  while(nchar(left)>0) {
    if(convertWidth(grobWidth(textGrob(left,gp=gp)),'mm',valueOnly=TRUE)<=width){ans<-c(ans,left);break}
    lens <- vapply(seq_len(nchar(left)),function(k)convertWidth(grobWidth(textGrob(substr(left,1,k),gp=gp)),'mm',valueOnly=TRUE),numeric(1))
    end <- max(which(lens<=width));cut <- gregexpr('[- /]',substr(left,1,end))[[1]];cut <- cut[cut>floor(end*.55)]
    if(length(cut))end <- max(cut)
    ans <- c(ans,trimws(substr(left,1,end)));left <- trimws(substr(left,end+1,nchar(left)))
  }
  stopifnot(length(ans)<=2);ans
}
draw <- function(record=FALSE) {
  grid.newpage();pushViewport(viewport(width=unit(width_mm,'mm'),height=unit(height_mm,'mm'),xscale=c(0,width_mm),yscale=c(height_mm,0)))
  txt <- function(s,x,y,size=7,face='plain',col='#263541',just='left')grid.text(s,x=unit(x,'native'),y=unit(y,'native'),just=just,gp=gpar(fontfamily='Helvetica',fontsize=size,fontface=face,col=col))
  for(j in 1:3)txt(titles[j],xs[j]+hw/2,7,8,'bold',just='centre')
  cursor <- 12;panels <- list();labels <- list();gidx <- 0
  for(g in groups) {
    idx <- which(r$group==g);nr <- length(idx)
    # A category heading is a structural label; scientific details stay in Handoff.
    txt(g,5,cursor+3,8,'bold');top <- cursor+7
    for(j in 1:3) {
      mat <- do.call(rbind,lapply(r$common_id[idx],function(id){a<-d[d$common_id==id & d$donor==donors[j],];a<-a[order(a$t),];stopifnot(nrow(a)>1);a$display_z}))
      grid.raster(matrix(colors(mat),nrow=nr),x=unit(xs[j]+hw/2,'native'),y=unit(top+nr*rowh/2,'native'),width=unit(hw,'mm'),height=unit(nr*rowh,'mm'),interpolate=FALSE)
      pushViewport(viewport(x=unit(xs[j]+hw/2,'native'),y=unit(top+nr*rowh/2,'native'),width=unit(hw,'mm'),height=unit(nr*rowh,'mm')))
      lo <- deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE);hi <- deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE);popViewport()
      panels[[length(panels)+1]] <- list(id=paste0('group',gidx,'donor',j),bbox_pt=c(lo$x,lo$y,hi$x,hi$y)*72,grid_id='heatmap',row_start=gidx,row_stop=gidx+1,col_start=j-1,col_stop=j)
      for(k in seq_along(idx)) {
        ii <- idx[k];q <- r[ii,paste0(donors[j],'_q')];delta <- r[ii,paste0(donors[j],'_end_difference')]
        txt(paste0(if(delta>0)'U' else if(delta<0)'D' else '0',if(q<settings$dynamic_alpha)'*' else ''),xs[j]+hw+3.7,top+(k-.5)*rowh,7,col=if(delta>0)'#B2182B' else '#2166AC',just='centre')
      }
    }
    for(k in seq_along(idx)) {
      ii <- idx[k];lab <- wrap_name(r$label[ii]);y <- top+(k-.5)*rowh
      for(l in seq_along(lab))txt(lab[l],5,y+(l-(length(lab)+1)/2)*2.65,7,if(r$is_SM[ii])'bold' else 'plain')
      labels[[length(labels)+1]] <- list(common_id=r$common_id[ii],full_label=r$label[ii],lines=lab)
    }
    cursor <- top+nr*rowh+3;gidx <- gidx+1
  }
  for(j in 1:3)for(t in c(0,.5,1))txt(format(t,trim=TRUE),xs[j]+hw*t,cursor+1,7,just='centre')
  txt('Relative pseudotime',177,cursor+7,8,just='centre')
  txt('Row Z',169,cursor+18,7,just='right')
  grid.raster(matrix(pal,nrow=1),x=unit(192,'native'),y=unit(cursor+18,'native'),width=unit(37,'mm'),height=unit(2.5,'mm'),interpolate=FALSE)
  for(z in c(-2,0,2))txt(as.character(z),173.5+(z+2)/4*37,cursor+23,7,just='centre')
  if(record) {
    write_json(list(schema_version=1L,backend='R grid deviceLoc',figure=list(width_pt=width_mm/25.4*72,height_pt=height_mm/25.4*72),panels=panels),file.path(qa,'figure4e_layout.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
    write_json(labels,file.path(qa,'label_layout.json'),auto_unbox=TRUE,pretty=TRUE)
  }
  popViewport()
}
dir.create(file.path(panel_dir,'FinalPdf'),showWarnings=FALSE)
pdf(NULL,width=width_mm/25.4,height=height_mm/25.4,family='Helvetica',useDingbats=FALSE);draw(TRUE);dev.off()
pdf(file.path(panel_dir,'FinalPdf/figure4e.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Helvetica',useDingbats=FALSE);draw();dev.off()
ragg::agg_png(file.path(qa,'figure4e_preview.png'),width=width_mm,height=height_mm,units='mm',res=160,background='white');draw();dev.off()
capture.output(sessionInfo(),file=file.path(qa,'R_session.txt'))
