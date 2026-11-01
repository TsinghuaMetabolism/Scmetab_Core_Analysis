# R rendering of collaborator box + reflected KDE method on current T3 t_graph.
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
source(file.path(dirname(normalizePath(script)),'plot_common.R'))
settings <- fromJSON(args[2])

dd <- file.path(panel_dir,'OriginalData');rd <- function(n)fread(file.path(dd,n))
c <- rd('cells.tsv');b <- rd('group_summary.tsv');v <- rd('distribution_curves.tsv');bins <- rd('distribution_bins.tsv')
s <- fromJSON(file.path(dd,'statistics.json'))
# Independent rank test and type7 quartiles; all projected cells participate.
test <- wilcox.test(c[group=='Myr']$t_graph_Myr,c[group=='Ctr']$t_graph_Myr,exact=FALSE,correct=TRUE)
stopifnot(abs(unname(test$statistic)-s$U_Myr)<1e-8,abs(test$p.value-s$P)<max(1e-15,s$P*1e-10))
for(g in c('Ctr','Myr')) {
  raw <- c[group==g]$t_graph_Myr;row <- b[group==g]
  stopifnot(max(abs(quantile(raw,c(.25,.5,.75),type=7)-as.numeric(row[,.(q1,median,q3)])))<1e-12)
  z <- v[group==g];bw <- s$bandwidth
  den <- vapply(z$t,function(t)mean(dnorm((t-raw)/bw)+dnorm((t+raw)/bw)+dnorm((t-2+raw)/bw))/bw,numeric(1))
  den <- den/sum(diff(z$t)*(head(den,-1)+tail(den,-1))/2)
  stopifnot(max(abs(den-z$density))<1e-12,abs(sum(bins[group==g]$percent)-100)<1e-10)
}
write_json(list(status='PASS',independent_backend='R wilcox.test + quantile(type=7) + reflected Gaussian KDE',n=as.list(s$n)),file.path(dd,'independent_statistics_verification.json'),auto_unbox=TRUE,pretty=TRUE)
colours <- c(Ctr='#2878B5',Myr='#AC577C')
b[,y:=ifelse(group=='Ctr',1,0)]
top <- ggplot(b)+
  geom_segment(aes(x=whisker_low,xend=whisker_high,y=y,yend=y),linewidth=.4,colour='#505962')+
  geom_segment(aes(x=whisker_low,xend=whisker_low,y=y-.12,yend=y+.12),linewidth=.4,colour='#505962')+
  geom_segment(aes(x=whisker_high,xend=whisker_high,y=y-.12,yend=y+.12),linewidth=.4,colour='#505962')+
  geom_rect(aes(xmin=q1,xmax=q3,ymin=y-.22,ymax=y+.22,fill=group,colour=group),alpha=.55,linewidth=.4)+
  geom_segment(aes(x=median,xend=median,y=y-.22,yend=y+.22),linewidth=.7,colour='#263341')+
  annotate('segment',x=1.025,xend=1.04,y=0,yend=0,linewidth=.4)+
  annotate('segment',x=1.04,xend=1.04,y=0,yend=1,linewidth=.4)+
  annotate('segment',x=1.025,xend=1.04,y=1,yend=1,linewidth=.4)+
  annotate('text',x=1.062,y=.5,label=s$stars,hjust=0,size=2.8,family='Helvetica')+
  scale_fill_manual(values=colours)+scale_colour_manual(values=colours)+
  scale_y_continuous(breaks=c(0,1),labels=c('Myr','Ctr'),limits=c(-.55,1.55),expand=expansion(mult=0))+
  scale_x_continuous(breaks=c(0,.25,.5,.75,1))+
  coord_cartesian(xlim=c(0,1),expand=FALSE,clip='off')+labs(x=NULL,y=NULL)+base_theme()+
  theme(legend.position='none',axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),plot.margin=margin(2,10,0,2,'mm'))
step <- rbindlist(lapply(c('Ctr','Myr'),function(g) {
  a <- bins[group==g];data.table(group=g,t=c(a$left,1),percent=c(a$percent,tail(a$percent,1)))
}))
bottom <- ggplot(v,aes(t,percent_per_bin_width,colour=group,fill=group))+
  geom_step(data=step,aes(t,percent,colour=group),inherit.aes=FALSE,linewidth=.25,alpha=.4,direction='hv')+
  geom_area(alpha=.12,position='identity',colour=NA)+geom_line(linewidth=.65)+
  scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
  scale_x_continuous(breaks=c(0,.25,.5,.75,1),limits=c(0,1),expand=expansion(mult=0))+
  scale_y_continuous(expand=expansion(mult=c(0,.06)))+
  labs(x='Pseudotime',y='Cells (% per 0.05)')+base_theme()+
  theme(legend.position='none',axis.text.x=element_text(margin=margin(t=5,unit='pt')),
        plot.margin=margin(2,10,2,2,'mm'))
render_grid(list(top,bottom),panel_dir,'figure5f',2,1,112,112,row_weights=c(1,2.3),legend=colours)
