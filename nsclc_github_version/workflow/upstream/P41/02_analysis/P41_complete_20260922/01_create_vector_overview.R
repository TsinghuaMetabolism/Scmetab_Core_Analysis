# Vector-only overview: reuse the local CSV plotting scripts and their R grobs.
# No PNG embedding, no change to panel data, no overwrite of panel figures.
PAGE_WIDTH_IN<-32;PAGE_HEIGHT_IN<-44
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root<-dirname(normalizePath(script))
suppressPackageStartupMessages({library(grid);library(ggplot2);library(patchwork)})
index<-read.csv(file.path(root,'panel_script_index.csv'))
plots<-list()
# Scratch graphics device for pheatmap layout, never creates an output file.
grDevices::pdf(NULL,width=15,height=14)
for(i in seq_len(nrow(index))){
 panel<-index$panel[i];paneldir<-file.path(root,panel)
 env<-new.env(parent=globalenv());env$captured<-list()
 env$ggsave<-eval(quote(function(filename,plot,...){
   if(grepl('[.]pdf$',filename))captured[[basename(filename)]]<<-plot
   invisible(NULL)
 }),env)
 env$write.csv<-function(...)invisible(NULL)
 env$cairo_pdf<-function(...)invisible(NULL)
 env$agg_png<-function(...)invisible(NULL)
 env$dev.off<-function(...)invisible(NULL)
 env$print<-function(...)invisible(NULL)
 code<-readLines(file.path(paneldir,index$R_script[i]))
 # Standalone scripts normally use --file to locate themselves; supply panel path here.
 code[grep('source_file<-',code)]<-paste0('setwd(',encodeString(paneldir,quote='"'),')')
 eval(parse(text=code),env)
 if(panel%in%c('panel_c','panel_e')){
   plots[[panel]]<-env$draw()$gtable
 }else{
   for(n in names(env$captured)){
     p<-env$captured[[n]]
     name<-if(panel=='panel_g')if(grepl('g1',n))'panel_g1' else 'panel_g2' else panel
     plots[[name]]<-if(inherits(p,'patchwork'))patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)
   }
 }
}
grDevices::dev.off();setwd(root)
order<-c('panel_a','panel_b','panel_c','panel_d','panel_e','panel_g1','panel_g2','panel_h')
titles<-c('a  |  Batch integration','b  |  Osi reference trajectory','c  |  Dynamic heatmap: K = 4','d  |  KEGG enrichment','e  |  MWU + effect size: selected metabolites','g1  |  Fixed-reference projection','g2  |  Reference interval composition','h  |  Significant choline-pathway dynamics')
stopifnot(all(order%in%names(plots)))
grDevices::cairo_pdf(file.path(root,'overview.pdf'),width=PAGE_WIDTH_IN,height=PAGE_HEIGHT_IN,onefile=TRUE)
grid.newpage();grid.rect(gp=gpar(fill='#F0F4F8',col=NA))
grid.text('P41 NSCLC PDO | Eight-panel overview',x=.035,y=.978,just='left',gp=gpar(fontsize=32,fontface='bold',col='#17324D'))
grid.text('Vector PDF | CSV inputs and named R scripts supplied',x=.035,y=.957,just='left',gp=gpar(fontsize=17,col='#52687C'))
for(i in seq_along(order)){
 column<-(i-1)%%2;row<-(i-1)%/%2
 pushViewport(viewport(x=.025+column*.49+.235,y=.916-row*.224-.099,width=.47,height=.213))
 grid.roundrect(r=unit(.008,'npc'),gp=gpar(fill='white',col='#DFE6EC'))
 grid.text(titles[i],x=.025,y=.967,just='left',gp=gpar(fontsize=18,fontface='bold',col='#17324D'))
 pushViewport(viewport(x=.5,y=.455,width=.95,height=.86))
 grid.draw(plots[[order[i]]]);popViewport(2)
}
grid.text('Pseudobulk display: same contrast BH-FDR < 0.05 and |rank-biserial| >= 0.30.',x=.035,y=.02,just='left',gp=gpar(fontsize=15,col='#52687C'))
grDevices::dev.off()
message('Created vector overview.pdf; all original panel inputs and outputs unchanged.')
