# Portable R renderer. The geometry measurement follows the existing grid renderer.
suppressPackageStartupMessages({library(ggplot2);library(data.table);library(jsonlite);library(grid)})
day_colors <- c(Day0='#2878B5', Day3='#E69F00', Day9='#8064A2')
group_colors <- c(C0='#777777', M0='#17A6A6')
region_colors <- c(Early='#D7E7F2', Middle='#84AFCB', Late='#2E5978')
base_theme <- function(square=FALSE) {
  theme_classic(base_size=7, base_family='Helvetica') +
    theme(axis.title=element_text(size=7), axis.text=element_text(size=6),
          axis.line=element_line(linewidth=.25), axis.ticks=element_line(linewidth=.25),
          plot.title=element_text(size=7.5, face='bold', margin=margin(b=4)),
          plot.margin=margin(2,2,2,2,'mm'), legend.position='bottom',
          legend.text=element_text(size=6), legend.title=element_text(size=6),
          legend.key.size=unit(2.6,'mm'), legend.spacing.x=unit(1,'mm'),
          legend.background=element_blank(),legend.key=element_blank(),
          aspect.ratio=if(square) 1 else NULL)
}

root_star <- function(x,y,xrange,yrange) {
  # Ten alternating vertices form a filled five-point star. Its centre is
  # exactly the frozen root; relative radii preserve its shape on square axes.
  angle <- pi/2 + (0:9)*pi/5
  radius <- rep(c(1,(3-sqrt(5))/2),5)
  vertices <- data.frame(x=x+cos(angle)*radius*diff(xrange)*.024,
                         y=y+sin(angle)*radius*diff(yrange)*.024)
  geom_polygon(data=vertices,aes(x,y),inherit.aes=FALSE,
               fill='black',colour='black',linewidth=.12,show.legend=FALSE)
}

render_grid <- function(plots, panel_dir, name, nr, nc, width_mm, height_mm,
                        header='', columns=NULL, footer='', square=FALSE, legend=NULL, fixed_panel_mm=NULL,
                        row_weights=rep(1,nr)) {
  stopifnot(length(plots)==nr*nc)
  out <- file.path(panel_dir, 'FinalPdf')
  qa <- file.path(panel_dir, 'OriginalData/qa')
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  dir.create(qa, recursive=TRUE, showWarnings=FALSE)
  # Open the measurement device before ggplotGrob queries font metrics, so R
  # never opens its implicit Rplots.pdf device in the caller's directory.
  pdf(NULL,width=width_mm/25.4,height=height_mm/25.4,family='Helvetica',useDingbats=FALSE)
  grobs <- lapply(plots, ggplotGrob)
  if(!is.null(fixed_panel_mm)) grobs <- lapply(grobs,function(g) {
    k <- which(g$layout$name=='panel');stopifnot(length(k)==1)
    g$widths[g$layout$l[k]] <- unit(fixed_panel_mm,'mm')
    g$heights[g$layout$t[k]] <- unit(fixed_panel_mm,'mm')
    g
  })
  widths <- do.call(unit.pmax, lapply(grobs, function(g) g$widths))
  heights <- do.call(unit.pmax, lapply(grobs, function(g) g$heights))
  grobs <- lapply(grobs, function(g) {g$widths<-widths;g$heights<-heights;g})
  top <- if(nzchar(header)) 9 else 0
  column_height <- if(is.null(columns)) 0 else 6
  bottom <- if(!is.null(legend)) 6 else if(nzchar(footer)) 10 else 0
  start <- 3
  stopifnot(length(row_weights)==nr,all(row_weights>0))
  row_heights <- unit.c(unit(c(top,column_height),'mm'),unit(row_weights,'null'),unit(bottom,'mm'))
  draw <- function(measure=FALSE) {
    grid.newpage()
    pushViewport(viewport(name='figure-page', layout=grid.layout(length(row_heights), nc, heights=row_heights)))
    if(nzchar(header)) {
      pushViewport(viewport(layout.pos.row=1, layout.pos.col=1:nc))
      grid.text(header, x=unit(3,'mm'), just='left', gp=gpar(fontfamily='Helvetica', fontsize=9, fontface='bold'))
      popViewport()
    }
    if(!is.null(columns)) for(j in seq_len(nc)) {
      pushViewport(viewport(layout.pos.row=2, layout.pos.col=j))
      grid.text(columns[j], gp=gpar(fontfamily='Helvetica', fontsize=8, fontface='bold'))
      popViewport()
    }
    panels <- list()
    for(j in seq_along(grobs)) {
      pushViewport(viewport(name=paste0('slot',j), layout.pos.row=start+(j-1)%/%nc, layout.pos.col=1+(j-1)%%nc))
      grid.draw(grobs[[j]])
      if(measure) {
        grid.force()
        v <- grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)
        idx <- which(grepl('^panel[.]',v$name) & grepl(paste0('slot',j,'::'),paste0(v$vpPath,'::'),fixed=TRUE))
        stopifnot(length(idx)==1)
        path <- strsplit(paste(v$vpPath[idx],v$name[idx],sep='::'),'::',fixed=TRUE)[[1]]
        upViewport(0);downViewport(do.call(vpPath,as.list(path[-1])))
        lo <- deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE)
        hi <- deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
        bb <- c(lo$x,lo$y,hi$x,hi$y)*72
        if(square) stopifnot(abs((bb[3]-bb[1])-(bb[4]-bb[2])) < 1e-6)
        panels[[j]] <- list(id=paste0(name,'_',j), bbox_pt=bb, grid_id='main',
                            row_start=(j-1)%/%nc, row_stop=1+(j-1)%/%nc,
                            col_start=(j-1)%%nc, col_stop=1+(j-1)%%nc)
        upViewport(0);downViewport('figure-page')
      } else popViewport()
    }
    pushViewport(viewport(layout.pos.row=length(row_heights),layout.pos.col=1:nc))
    if(nzchar(footer)) grid.text(footer, x=unit(3,'mm'), just='left', gp=gpar(fontfamily='Helvetica', fontsize=6, lineheight=1.25))
    if(!is.null(legend)) for(j in seq_along(legend)) {
      x <- unit(.5,'npc')+unit((j-(length(legend)+1)/2)*16-4,'mm')
      grid.rect(x=x,y=.5,width=unit(2.8,'mm'),height=unit(2.8,'mm'),gp=gpar(fill=legend[j],col=NA))
      grid.text(names(legend)[j],x=x+unit(2.5,'mm'),just='left',gp=gpar(fontfamily='Helvetica',fontsize=7))
    }
    upViewport(0)
    if(measure) {
      manifest <- list(schema_version=1L, backend='R grid deviceLoc',
                       figure=list(width_pt=width_mm/25.4*72,height_pt=height_mm/25.4*72),
                       panels=panels, square_plot_areas=square)
      manifest_file <- file.path(qa,paste0(name,'_layout.json'))
      write_json(manifest,manifest_file,auto_unbox=TRUE,pretty=TRUE,digits=16)

    }
  }
  # Measure and enforce alignment before opening the final output file.
  draw(TRUE);dev.off()
  pdf(file.path(out,paste0(name,'.pdf')),width=width_mm/25.4,height=height_mm/25.4,family='Helvetica',useDingbats=FALSE)
  draw();dev.off()
  ragg::agg_png(file.path(qa,paste0(name,'_preview.png')),width=width_mm,height=height_mm,units='mm',res=160,background='white')
  draw();dev.off()
  capture.output(sessionInfo(), file=file.path(qa,'R_session.txt'))
}
