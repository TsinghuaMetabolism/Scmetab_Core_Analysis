# Redraw saved complete reference coordinates, with no fitting or displacement.
args <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
workflow <- dirname(dirname(dirname(normalizePath(script))))
source(file.path(workflow,'scripts/plot_common.R'))
settings <- fromJSON(args[2])

name <- basename(panel_dir)
rd <- function(n) fread(file.path(panel_dir,'OriginalData',n))
cells <- rd('cells.tsv');nodes <- rd('centers.tsv');edges <- rd('center_edges.tsv')
m <- fromJSON(file.path(panel_dir,'OriginalData/model.json'))
stopifnot(nrow(edges)==nrow(nodes)-1)
cells[,day:=factor(round,levels=c(0,1,3),labels=c('Day0','Day3','Day9'))]
segments <- data.table(x=nodes$DDRTree1[match(edges$from,nodes$center_id)],
                       y=nodes$DDRTree2[match(edges$from,nodes$center_id)],
                       xend=nodes$DDRTree1[match(edges$to,nodes$center_id)],
                       yend=nodes$DDRTree2[match(edges$to,nodes$center_id)])
root <- cells[cell_id==m$root_cell];stopifnot(nrow(root)==1)
base <- ggplot() + geom_segment(data=segments,aes(x,y,xend=xend,yend=yend),linewidth=.45,colour='#263341') +
  labs(x='DDRTree 1',y='DDRTree 2') + base_theme(TRUE)
points <- function(aes) geom_point(data=cells,mapping=aes,shape=21,size=1.05,stroke=.12,alpha=.85)
fixed_mm <- NULL
xl <- range(c(cells$DDRTree1,nodes$DDRTree1));yl <- range(c(cells$DDRTree2,nodes$DDRTree2))
if(name=='figure4c') {
  view <- fromJSON(file.path(panel_dir,'OriginalData/reference_display.json'))
  xl <- view$xlim;yl <- view$ylim;fixed_mm <- view$plot_area_mm
  stopifnot(max(abs(c(root$DDRTree1,root$DDRTree2)-view$root_xy))<1e-13)
  base <- base+coord_cartesian(xlim=xl,ylim=yl,expand=FALSE)
}
root_layer <- root_star(root$DDRTree1,root$DDRTree2,xl,yl)
p1 <- base + points(aes(DDRTree1,DDRTree2,fill=day)) + root_layer +
  scale_fill_manual(values=day_colors,name=NULL,labels=c('Day 0','Day 3','Day 9'))
p2 <- base + points(aes(DDRTree1,DDRTree2,fill=Pseudotime01)) + root_layer +
  scale_fill_gradientn(colours=c('#440154','#31688E','#35B779','#FDE725'),limits=c(0,1),breaks=c(0,.5,1),name='Pseudotime',
                       guide=guide_colourbar(barwidth=unit(28,'mm'),barheight=unit(2,'mm'),title.position='top',title.hjust=.5)) +
  theme(legend.text=element_text(size=6,margin=margin(t=3,unit='pt')))
render_grid(list(p1,p2),panel_dir,name,1,2,183,96,square=TRUE,fixed_panel_mm=fixed_mm)
