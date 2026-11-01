suppressPackageStartupMessages(library(data.table))
a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-dirname(normalizePath(sub('^--file=','',a[1])))
r<-function(p)fread(file.path(root,p))
input<-r('input_validation.csv')
a_cells<-r('panel_a/tables/plot_data.csv')
b_cells<-r('panel_b/tables/plot_cells.csv')
b_screen<-r('panel_b/inputs/ordering_feature_screen.csv')
c_stats<-r('panel_c/tables/all_features_GAM_significance.csv')
c_heat<-r('panel_c/tables/heatmap_rowZ.csv')
c_clusters<-r('panel_c/tables/cluster_counts.csv')
d_query<-r('panel_d/tables/GAM_significant_annotated_query_features.csv')
d_plot<-r('panel_d/tables/plot_data.csv')
e_summary<-r('panel_e/tables/analysis_summary.csv')
heldout<-r('panel_a/tables/heldout_leakage_checks.csv')
comparison<-r('tables/previous_vs_endpoint_summary.csv')
expected_c<-merge(r('panel_c/inputs/features.csv')[!is.na(annotation)&nzchar(trimws(annotation)),.(feature_id)],
                  c_stats[association_BH_FDR_Osi<.05,.(feature_id)],by='feature_id')
checks<-data.table(
 check=c('input checks passed','panel A has 1068 labeled cells',
         'panel B has 1068 cells','panel B uses 159 ordering features',
         'panel B has one natural state','panel C tests 3338 features',
         'panel C heatmap equals significant annotated set','panel C has four clusters',
         'panel D query equals panel C heatmap','panel D sorted by BH FDR',
         'panel E BH family is 3338 x 4','panel E MWU formula matches wilcox',
         'heldout split leakage checks passed','old-vs-new comparison exported',
         'five main PDFs exist','five main PNGs exist'),
 passed=c(all(input$passed),nrow(a_cells)==1068L&all(nzchar(a_cells$cell_type)),
          nrow(b_cells)==1068L,sum(b_screen$selected)==159L,
          uniqueN(b_cells$Monocle_state)==1L,nrow(c_stats)==3338L,
          setequal(c_heat$feature_id,expected_c$feature_id),nrow(c_clusters)==4L,
          setequal(d_query$feature_id,c_heat$feature_id),
          nrow(d_plot)<2L||all(diff(d_plot$p.adjust)>=0),
          e_summary[metric=='BH_family_size',as.numeric(value)]==13352,
          e_summary[metric=='MWU_formula_vs_wilcox_max_difference',as.numeric(value)]<1e-10,
          all(heldout$passed),nrow(comparison)==11L,
          all(file.exists(file.path(root,c('panel_a/figures/panel_a_UMAP_before_after.pdf',
            'panel_b/figures/panel_b_Osi_reference_trajectory.pdf',
            'panel_c/figures/panel_c_K4_dynamic_heatmap.pdf',
            'panel_d/figures/panel_d_KEGG_enrichment.pdf',
            'panel_e/figures/panel_e_pseudobulk_heatmap.pdf')))),
          all(file.exists(file.path(root,c('panel_a/figures/panel_a_UMAP_before_after.png',
            'panel_b/figures/panel_b_Osi_reference_trajectory.png',
            'panel_c/figures/panel_c_K4_dynamic_heatmap.png',
            'panel_d/figures/panel_d_KEGG_enrichment.png',
            'panel_e/figures/panel_e_pseudobulk_heatmap.png')))))
)
fwrite(checks,file.path(root,'validation.csv'))
print(checks);stopifnot(all(checks$passed))
cat('All',nrow(checks),'validation checks passed.\n')
