# Rebuild eight overview figures plus the treatment-split UMAP supplement.
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-dirname(normalizePath(script))
panels<-paste0('panel_',c('a','b','c','d','e','g','h'))
script_index<-read.csv(file.path(root,'panel_script_index.csv'))
if(!'--preview-only'%in%commandArgs(TRUE))for(p in panels){
 status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(root,p,script_index$R_script[match(p,script_index$panel)])))
 if(status!=0)stop('Failed: ',p)
}
if(!'--preview-only'%in%commandArgs(TRUE)){
 status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(root,'panel_a/scripts/02_plot_integrated_UMAP_split_by_treatment.R')))
 if(status!=0)stop('Failed: panel_a treatment-split UMAP')
}
suppressPackageStartupMessages({library(grid);library(png)})
files<-c('panel_a/figures/panel_a_UMAP_before_after.png','panel_b/figures/panel_b_Osi_reference_trajectory.png','panel_c/figures/panel_c_K4_dynamic_heatmap.png','panel_d/figures/panel_d_KEGG_enrichment.png','panel_e/figures/panel_e_pseudobulk_heatmap.png','panel_g/figures/panel_g1_fixed_reference_projection.png','panel_g/figures/panel_g2_state_composition.png','panel_h/figures/panel_h_significant_choline_pathway_dynamics.png')
titles<-c('a · Batch integration','b · Osi reference trajectory','c · Dynamic heatmap | K = 4','d · KEGG pathway enrichment','e · MWU + effect size | selected features','g1 · Fixed-reference projection','g2 · Reference interval composition','h · Significant choline-pathway dynamics')
stopifnot(all(file.exists(file.path(root,files))))
png(file.path(root,'overview.png'),width=2400,height=3200,res=150,type='cairo')
grid.newpage();grid.rect(gp=gpar(fill='#EDF2F7',col=NA));grid.text('P41 NSCLC PDO | Eight-figure edition',x=.045,y=.975,just='left',gp=gpar(fontsize=24,fontface='bold',col='#17324D'))
for(i in seq_along(files)){
 co<-(i-1)%%2;ro<-(i-1)%/%2
 pushViewport(viewport(x=.025+co*.49+.235,y=.91-ro*.224-.1,width=.47,height=.21))
 grid.roundrect(r=unit(.015,'npc'),gp=gpar(fill='white',col='#DFE6EC'))
 grid.text(titles[i],x=.035,y=.935,just='left',gp=gpar(fontsize=11,fontface='bold',col='#17324D'))
 im<-readPNG(file.path(root,files[i]));ar<-dim(im)[2]/dim(im)[1];boxar<-(2400*.47)/(3200*.21);w<-.93;h<-w*boxar/ar
 if(h>.79){h<-.79;w<-h/boxar*ar}
 grid.raster(im,x=.5,y=.46,width=w,height=h,interpolate=TRUE);popViewport()
}
grid.text('CSV inputs | Named R scripts | Processed tables | PDF and PNG',x=.045,y=.025,just='left',gp=gpar(fontsize=12,col='#52687C'));dev.off()
# Small data/figure validation table, no extra diagnostic figures.
cdata<-read.csv(file.path(root,'panel_c/tables/heatmap_rowZ.csv'))
es<-read.csv(file.path(root,'panel_e/tables/analysis_summary.csv'))
dd<-read.csv(file.path(root,'panel_d/tables/plot_data.csv'))
ss<-read.csv(file.path(root,'panel_h/inputs/all_features_GAM_significance.csv'))
aa<-read.csv(file.path(root,'panel_h/inputs/annotated_feature_KEGG_pathway_mapping.csv'))
expected<-intersect(aa$feature_id,ss$feature_id[ss$association_BH_FDR_Osi<.05])
checks<-data.frame(check=c('Eight main PDFs plus split and time-resolved UMAP supplements','Eight main PNGs plus two UMAP supplements and overview','Panel C four clusters','Panel C feature set equals shared GAM significant annotations','Panel C FDR equals panel H','Panel D increasing FDR top to bottom','Panel E all feature-by-day BH tests','Panel E MWU independent formula check'),passed=c(sum(basename(list.files(root,pattern='[.]pdf$',recursive=TRUE))!='overview.pdf')==10,length(list.files(root,pattern='[.]png$',recursive=TRUE))==11,length(unique(cdata$cluster))==4,setequal(cdata$feature_id,expected)&&!anyDuplicated(cdata$feature_id),max(abs(cdata$q_value-ss$association_BH_FDR_Osi[match(cdata$feature_id,ss$feature_id)]))<1e-12,all(diff(dd$p.adjust)>=0),es$value[es$metric=='BH_family_size']==es$value[es$metric=='eligible_features']*4,es$value[es$metric=='MWU_formula_vs_wilcox_max_difference']<1e-10))
write.csv(checks,file.path(root,'validation.csv'),row.names=FALSE);stopifnot(all(checks$passed))
cat('Finished: eight overview figures and the split UMAP supplement; all checks passed.\n')
