# Native Ctr reference plus saved Myr projection; frozen numerical results retained.
args <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
workflow <- dirname(dirname(dirname(normalizePath(script))))
source(file.path(workflow,'scripts/plot_common.R'))
settings <- fromJSON(args[2])

rd <- function(n) fread(file.path(panel_dir,'OriginalData',n))
d <- rd('cells.tsv');nodes <- rd('centers.tsv');learned <- rd('center_edges.tsv')
reference <- rd('reference_cells.tsv');display <- rd('display_cells.tsv')
view <- fromJSON(file.path(panel_dir,'OriginalData/reference_display.json'))
stopifnot(all(d$group %in% c('Ctr','Myr')),all(table(d$group)>0))
stopifnot(nrow(learned)==nrow(nodes)-1,identical(d$cell_id,display$cell_id))
d[,`:=`(plot_x=display$plot_x,plot_y=display$plot_y)]
stopifnot(max(abs(as.matrix(d[group=='Myr',.(plot_x,plot_y)])-
                  as.matrix(d[group=='Myr',.(projection_x_Myr,projection_y_Myr)])))<1e-13)
native <- reference[match(d[group=='Ctr',cell_id],cell_id)]
stopifnot(max(abs(as.matrix(d[group=='Ctr',.(plot_x,plot_y)])-
                  as.matrix(native[,.(DDRTree1,DDRTree2)])))<1e-13)
edges <- data.table(x=nodes$DDRTree1[match(learned$from,nodes$center_id)],
                    y=nodes$DDRTree2[match(learned$from,nodes$center_id)],
                    xend=nodes$DDRTree1[match(learned$to,nodes$center_id)],
                    yend=nodes$DDRTree2[match(learned$to,nodes$center_id)])
d[,day:=factor(actual_day,levels=c(0,3,9),labels=c('Day0','Day3','Day9'))]
xl <- view$xlim;yl <- view$ylim
base <- ggplot() + geom_segment(data=edges,aes(x,y,xend=xend,yend=yend),linewidth=.45,colour='#263341') +
  coord_cartesian(xlim=xl,ylim=yl,expand=FALSE) + base_theme(TRUE) + labs(x='DDRTree 1',y='DDRTree 2')
gray <- function(x) geom_point(data=x,aes(plot_x,plot_y),colour='#A8A8A8',size=.55,alpha=.35)
r <- reference[cell_id==view$root_cell];stopifnot(nrow(r)==1)
stopifnot(max(abs(c(r$DDRTree1,r$DDRTree2)-view$root_xy))<1e-13)
root <- root_star(r$DDRTree1,r$DDRTree2,xl,yl)
p <- base + gray(d[group=='Ctr']) +
  geom_point(data=d[group=='Myr'],aes(plot_x,plot_y,fill=day),shape=21,size=1.0,alpha=.85,stroke=.1) +
  root + scale_fill_manual(values=day_colors,name=NULL,labels=c('Day 0','Day 3','Day 9')) + labs(title='Myr')
render_grid(list(p),panel_dir,'figure5e',1,1,91.5,96,square=TRUE,fixed_panel_mm=view$plot_area_mm)
