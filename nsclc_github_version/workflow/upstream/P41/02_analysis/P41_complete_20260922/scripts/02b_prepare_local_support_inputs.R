# Freeze P41 reference graph and existing Ridge prediction as CSV for standalone mapping.
suppressPackageStartupMessages({library(data.table);library(monocle);library(igraph)})
a<-grep('^--file=',commandArgs(),value=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',a[1]))))
d<-readRDS(file.path(root,'audit/analyzed.rds'));m<-d$metadata
# On a rerun, use preserved raw Ridge fields when available.
if('original_Ridge_1'%in%names(m)){m[,trajectory_1:=ifelse(batch=='B3',original_Ridge_1,trajectory_1)];m[,trajectory_2:=ifelse(batch=='B3',original_Ridge_2,trajectory_2)];m[,pseudotime:=ifelse(batch=='B3',original_Ridge_pseudotime,pseudotime)]}
m[,cohort:=ifelse(batch=='B3','Osi_query','Osi_reference')]
fwrite(m,file.path(root,'panel_g/inputs/original_Ridge_projection.csv'))
cds<-readRDS(file.path(root,'audit/Osi_Monocle2_model.rds'));g<-minSpanningTree(cds);k<-t(reducedDimK(cds));stopifnot(nrow(k)==vcount(g),ncol(k)==2)
nodes<-data.table(node_id=V(g)$name,x=k[,1],y=k[,2],degree=degree(g));edges<-as.data.table(as_edgelist(g,names=TRUE));setnames(edges,c('from','to'));edges[,edge_id:=.I];edges[,`:=`(x0=nodes$x[match(from,nodes$node_id)],y0=nodes$y[match(from,nodes$node_id)],x1=nodes$x[match(to,nodes$node_id)],y1=nodes$y[match(to,nodes$node_id)])]
ref<-m[batch!='B3'];earliest<-ref[which.min(pseudotime)];root_node<-nodes$node_id[which.min((nodes$x-earliest$trajectory_1)^2+(nodes$y-earliest$trajectory_2)^2)];stopifnot(nodes[node_id==root_node,degree]==1)
fwrite(nodes,file.path(root,'panel_g/inputs/reference_graph_nodes.csv'));fwrite(edges,file.path(root,'panel_g/inputs/reference_graph_edges.csv'));fwrite(data.table(root_node=root_node,root_cell=earliest$cell_id,root_native_pseudotime=earliest$pseudotime),file.path(root,'panel_g/inputs/reference_graph_root.csv'))
cat('Graph root',root_node,' tips',sum(nodes$degree==1),' branch nodes',sum(nodes$degree>2),'\n')
