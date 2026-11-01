# Install CSV-derived constrained mapping into the P41 analysis object; intensities stay identical.
suppressPackageStartupMessages(library(data.table))
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))));p<-file.path(root,'audit/analyzed.rds');d<-readRDS(p);old<-d$metadata;new<-fread(file.path(root,'panel_g/inputs/reference_and_projected_cells.csv'))
stopifnot(identical(old$cell_id,new$cell_id),max(abs(old$trajectory_1[old$batch!='B3']-new$trajectory_1[new$batch!='B3']))<1e-12,max(abs(old$pseudotime[old$batch!='B3']-new$pseudotime[new$batch!='B3']))<1e-12)
d$metadata<-new;saveRDS(d,p,compress=FALSE)
# Reference-only panel b keeps native graph values but adopts frozen K3 descriptive intervals.
fwrite(new[batch!='B3'],file.path(root,'panel_b/inputs/Osi_cells_and_trajectory.csv'))
cat('Installed',nrow(new),'cells; changed query coordinates',sum(new$constraint_active),'\n')
