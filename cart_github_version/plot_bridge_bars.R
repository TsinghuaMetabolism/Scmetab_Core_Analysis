# Same-run Bridge C0/M0 means and SEM.
args <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
source(file.path(dirname(normalizePath(script)),'plot_common.R'))
settings <- fromJSON(args[2])

stopifnot(basename(panel_dir)=='figure5d')
d <- fread(file.path(panel_dir,'OriginalData/group_summary.tsv'))
tests <- fread(file.path(panel_dir,'OriginalData/comparisons.tsv'))
d[,sem:=sd/sqrt(n_cells)]
d[,x:=match(donor,settings$bridge_donors)+ifelse(comparison_group=='C0',-.17,.17)]
stopifnot(all(d$n_cells>1))
stars <- function(q) ifelse(q<.0001,'****',ifelse(q<.001,'***',ifelse(q<.01,'**',ifelse(q<.05,'*','ns'))))
plots <- list();axes <- list()
for(target in settings$target_ids) {
  v <- d[target_id==target]
  z <- copy(tests[target_id==target])
  low <- min(v$mean-v$sem);high <- max(v$mean+v$sem)
  span <- max(high-low,.5)
  lower <- max(0,floor((low-.18*span)*2)/2)
  upper <- ceiling((high+.33*span)*2)/2
  full <- upper-lower
  v <- copy(v);v[,bar_base:=lower]
  z[,x:=match(donor,settings$bridge_donors)]
  z[,y:=high+.13*span]
  axes[[length(axes)+1]] <- data.table(target_id=target,label=v$label[1],lower=lower,upper=upper,
                                      truncated_axis=lower>0,error='sample SD / sqrt(n cells)')
  p <- ggplot(v) +
    geom_rect(aes(xmin=x-.145,xmax=x+.145,ymin=bar_base,ymax=mean,fill=comparison_group),colour=NA) +
    geom_errorbar(aes(x=x,ymin=mean-sem,ymax=mean+sem),width=.10,linewidth=.3,colour='#303030') +
    geom_segment(data=z,aes(x=x-.17,xend=x+.17,y=y,yend=y),linewidth=.25) +
    geom_segment(data=z,aes(x=x-.17,xend=x-.17,y=y-full*.025,yend=y),linewidth=.25) +
    geom_segment(data=z,aes(x=x+.17,xend=x+.17,y=y-full*.025,yend=y),linewidth=.25) +
    geom_text(data=z,aes(x=x,y=y+full*.07,label=stars(q_value)),size=2.5,family='Helvetica') +
    scale_fill_manual(values=group_colors,guide='none') +
    scale_x_continuous(breaks=1:3,labels=settings$bridge_donors,expand=c(0,0)) +
    scale_y_continuous(limits=c(lower,upper),breaks=pretty(c(lower,high),n=3),expand=c(0,0)) +
    labs(x=NULL,y=if(target==settings$target_ids[3]) 'TIC-normalized abundance' else NULL,title=v$label[1]) + base_theme() +
    theme(axis.ticks.x=element_blank(),axis.line.y=element_blank(),plot.title=element_text(size=7.5,face='bold'),
          plot.margin=margin(2,3,2,4,'mm'))
  if(target!=tail(settings$target_ids,1)) p <- p+theme(axis.text.x=element_blank())
  # Explicit breaks show the omitted zero-to-lower interval. All raw zeros
  # remain in both the means/SEM and the unchanged cell-level tests.
  if(lower>0) {
    cut <- lower+.045*full
    p <- p + geom_segment(aes(x=x-.155,xend=x+.155),y=cut-full*.015,yend=cut+full*.015,
                           colour='white',linewidth=1.15) +
      annotate('segment',x=.5,xend=.5,y=lower,yend=cut-full*.019,linewidth=.25) +
      annotate('segment',x=.5,xend=.5,y=cut+full*.021,yend=upper,linewidth=.25) +
      annotate('segment',x=.48,xend=.52,y=cut-full*.019,yend=cut+full*.005,linewidth=.25) +
      annotate('segment',x=.48,xend=.52,y=cut-full*.001,yend=cut+full*.023,linewidth=.25) +
      coord_cartesian(xlim=c(.5,3.5),clip='off')
  } else p <- p+annotate('segment',x=.5,xend=.5,y=lower,yend=upper,linewidth=.25)+coord_cartesian(xlim=c(.5,3.5))
  plots[[length(plots)+1]] <- p
}
fwrite(d,file.path(panel_dir,'OriginalData/bar_summary.tsv'),sep='\t')
fwrite(rbindlist(axes),file.path(panel_dir,'OriginalData/bar_axis_ranges.tsv'),sep='\t')
render_grid(plots,panel_dir,'figure5d',5,1,145,205,legend=group_colors)
