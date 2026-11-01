a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-dirname(normalizePath(sub('^--file=','',a[1])))
run<-function(p){
  cat('\n>>>',p,'\n')
  status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(root,p)))
  if(status!=0)stop('Failed: ',p)
}
run('00_prepare_inputs.R')
for(p in c('panel_a/scripts/01_plot_P39_UMAP_before_after_correction.R',
           'panel_a/scripts/02_plot_integrated_UMAP_split_by_batch.R',
           'panel_b/scripts/01_plot_Osi_DDRTree_group_and_pseudotime.R',
           'panel_b/scripts/02_plot_day_batch_and_stage_distribution.R',
           'panel_c/scripts/00_fit_single_condition_Gaussian_GAM.R',
           'panel_c/scripts/01_cluster_K4_and_plot_metabolite_dynamics.R'))run(p)
file.copy(file.path(root,'panel_c/tables/all_features_GAM_significance.csv'),
          file.path(root,'panel_d/inputs/all_features_GAM_significance.csv'),
          overwrite=TRUE)
run('panel_d/scripts/01_KEGG_ORA.R')
run('panel_e/scripts/01_test_timepoint_MWU_and_plot_pseudobulk.R')
run('03_compare_previous_analysis.R')
run('01_create_overview.R')
run('02_validate_results.R')
cat('\nP39 panels a-e analysis complete.\n')
