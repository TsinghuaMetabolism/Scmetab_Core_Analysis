suppressPackageStartupMessages({library(grid);library(png)})
a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-dirname(normalizePath(sub('^--file=','',a[1])))
files<-c('panel_a/figures/panel_a_UMAP_before_after.png',
         'panel_b/figures/panel_b_Osi_reference_trajectory.png',
         'panel_c/figures/panel_c_K4_dynamic_heatmap.png',
         'panel_d/figures/panel_d_KEGG_enrichment.png',
         'panel_e/figures/panel_e_pseudobulk_heatmap.png')
titles<-c('a · D0+D14 shared-anchor integration',
          'b · 159-feature D0-rooted trajectory',
          'c · Treated-cell Gaussian GAM | K = 4',
          'd · KEGG enrichment',
          'e · True-time MWU + effect size')
stopifnot(all(file.exists(file.path(root,files))))
draw<-function(){
 grid.newpage();grid.rect(gp=gpar(fill='#EDF2F7',col=NA))
 grid.text('P39 NSCLC PDO | Panels a–e',x=.04,y=.978,just='left',
           gp=gpar(fontsize=23,fontface='bold',col='#17324D'))
 grid.text('D0/D14 shared-anchor correction and endpoint-significant-feature trajectory',
           x=.04,y=.958,just='left',gp=gpar(fontsize=12,col='#52687C'))
 layout<-list(c(.27,.77,.44,.31),c(.74,.77,.44,.31),
              c(.27,.42,.44,.31),c(.74,.42,.44,.31),c(.50,.135,.70,.20))
 for(i in seq_along(files)){
   g<-layout[[i]];pushViewport(viewport(x=g[1],y=g[2],width=g[3],height=g[4]))
   grid.roundrect(r=unit(.012,'npc'),gp=gpar(fill='white',col='#DFE6EC'))
   grid.text(titles[i],x=.03,y=.95,just='left',
             gp=gpar(fontsize=10.5,fontface='bold',col='#17324D'))
   im<-readPNG(file.path(root,files[i]));ar<-dim(im)[2]/dim(im)[1]
   boxar<-(g[3]*2400)/(g[4]*3000);w<-.94;h<-w*boxar/ar
   if(h>.82){h<-.82;w<-h/boxar*ar}
   grid.raster(im,x=.5,y=.46,width=w,height=h,interpolate=TRUE);popViewport()
 }
 grid.text('Cell-level exploratory inference | one donor | D3=B1 only, D8=B2 only',
           x=.04,y=.012,just='left',gp=gpar(fontsize=11,col='#52687C'))
}
png(file.path(root,'overview.png'),width=2400,height=3000,res=150,type='cairo');draw();dev.off()
cairo_pdf(file.path(root,'overview.pdf'),width=16,height=20);draw();dev.off()
cat('Created overview.png and overview.pdf.\n')
