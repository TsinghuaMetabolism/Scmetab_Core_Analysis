# Five fixed SM rows by three donor columns; no data-dependent selection.
args <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
source(file.path(dirname(normalizePath(script)),'plot_common.R'))
settings <- fromJSON(args[2])

name <- basename(panel_dir)
stopifnot(name=='figure4g')
d <- fread(file.path(panel_dir,'OriginalData/abundance.tsv'))
s <- fread(file.path(panel_dir,'OriginalData/comparisons.tsv'))
donors <- settings$reference_donors
groups <- c('Day0','Day3','Day9')
targets <- settings$target_ids
qtext <- function(q) if(q < .001) sprintf('q = %.1e',q) else sprintf('q = %.3f',q)
plots <- list()
for(target in targets) {
  row_max <- max(d[target_id==target,abundance])
  step <- row_max*.15
  for(donor_id in donors) {
    v <- copy(d[donor==donor_id & target_id==target])
    v[,comparison_group:=factor(comparison_group,levels=groups)]
    z <- s[donor==donor_id & target_id==target]
    z <- z[c(1,3,2)] # short brackets first, spanning bracket last
    stopifnot(nrow(z)==3)
    top <- row_max + step*(nrow(z)+1.5)
    p <- ggplot(v,aes(comparison_group,abundance,fill=comparison_group)) +
      geom_violin(trim=TRUE,scale='width',bw='nrd0',adjust=1,width=.78,alpha=.7,linewidth=.25) +
      scale_fill_manual(values=day_colors,guide='none') +
      scale_x_discrete(labels=setNames(sub('Day','Day ',groups),groups)) +
      scale_y_continuous(limits=c(0,top),breaks=pretty(c(0,row_max),n=3),expand=expansion(mult=c(.02,.01))) +
      labs(x=NULL,y=if(donor_id==donors[1] && target==settings$target_ids[3]) 'TIC-normalized abundance' else NULL,
           title=if(donor_id==donors[1]) unique(v$label) else ' ') + base_theme()
    for(i in seq_len(nrow(z))) {
      x1 <- match(z$group_first[i],groups);x2 <- match(z$group_second[i],groups)
      y <- row_max + step*i
      p <- p + annotate('segment',x=x1,xend=x2,y=y,yend=y,linewidth=.25) +
        annotate('segment',x=c(x1,x2),xend=c(x1,x2),y=y-step*.18,yend=y,linewidth=.25) +
        annotate('text',x=(x1+x2)/2,y=y+step*.55,label=qtext(z$q_value[i]),size=1.95,family='Helvetica')
    }
    plots[[length(plots)+1]] <- p
  }
}
render_grid(plots,panel_dir,name,5,3,183,224,columns=settings$reference_labels)
