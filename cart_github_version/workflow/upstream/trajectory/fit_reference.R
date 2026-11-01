# CAR-T Ctr reference fitting and continuous-response spline basis.
# Supply study-specific settings separately; inputs are not included.

fit_reference <- function(x, metadata, parameters, baseline_stage, root_neighbors) {
  suppressPackageStartupMessages(library(monocle))
  x <- as.matrix(x) # ordering features x cells, already rank-normal transformed
  stopifnot(all(is.finite(x)), identical(colnames(x), rownames(metadata)),
            length(root_neighbors) == 3L, all(root_neighbors > 0),
            all(c('ncenter','lambda_mult','sigma','gamma','max_iter',
                  'expanded_max_iter','tolerance') %in% names(parameters)))
  pd <- new('AnnotatedDataFrame', data=metadata)
  fd <- new('AnnotatedDataFrame', data=data.frame(
    gene_short_name=rownames(x), row.names=rownames(x)))
  cds <- newCellDataSet(x, phenoData=pd, featureData=fd,
                       expressionFamily=VGAM::uninormal(), lowerDetectionLimit=0)
  cds <- setOrderingFilter(cds, rownames(x))
  normalized <- monocle:::normalize_expr_data(cds, norm_method='none', pseudo_expr=0)
  stopifnot(isTRUE(all.equal(unname(as.matrix(normalized)), unname(x), tolerance=0)))
  fit <- function(iterations) reduceDimension(cds, reduction_method='DDRTree',
      max_components=2, norm_method='none', pseudo_expr=0, relative_expr=FALSE,
      scaling=FALSE, auto_param_selection=FALSE,
      ncenter=parameters$ncenter, lambda=parameters$lambda_mult*ncol(x),
      sigma=parameters$sigma, param.gamma=parameters$gamma,
      maxIter=iterations, tol=parameters$tolerance)
  fitted <- fit(parameters$max_iter)
  objective <- as.numeric(fitted@auxOrderingData[['DDRTree']]$objective_vals)
  final_change <- function(v) tail(abs(diff(v))/pmax(abs(head(v,-1)), .Machine$double.eps),1)
  if (length(objective) >= parameters$max_iter && final_change(objective) >= parameters$tolerance) {
    fitted <- fit(parameters$expanded_max_iter)
    objective <- as.numeric(fitted@auxOrderingData[['DDRTree']]$objective_vals)
  }
  stopifnot(length(objective) > 1, all(is.finite(objective)))
  initial <- orderCells(fitted)
  states <- unique(as.character(pData(initial)$State))
  if (length(states) > 1) {
    candidates <- setNames(lapply(states, function(s) orderCells(initial, root_state=s)), states)
  } else {
    candidates <- list(forward=orderCells(initial, reverse=FALSE), reverse=orderCells(initial, reverse=TRUE))
  }
  candidates <- candidates[vapply(candidates, function(a)
    length(a@auxOrderingData[['DDRTree']]$root_cell)==1L, logical(1))]
  stopifnot(length(candidates) > 0)
  values <- t(x)
  roots <- do.call(rbind, lapply(names(candidates), function(key) {
    id <- candidates[[key]]@auxOrderingData[['DDRTree']]$root_cell
    distances <- rowSums(sweep(values,2,values[id,],'-')^2)
    stages <- metadata$round[order(distances)]
    data.frame(key=key, root_cell=id,
      early_primary=mean(head(stages,root_neighbors[1])==baseline_stage),
      early_second=mean(head(stages,root_neighbors[2])==baseline_stage),
      early_third=mean(head(stages,root_neighbors[3])==baseline_stage),
      mean_stage=mean(head(stages,root_neighbors[1])))
  }))
  roots <- roots[order(-roots$early_primary,-roots$early_second,-roots$early_third,roots$mean_stage),]
  roots$selected <- seq_len(nrow(roots))==1L
  rooted <- candidates[[roots$key[1]]]
  stopifnot(identical(reducedDimS(initial),reducedDimS(rooted)))
  list(cds=rooted, root_candidates=roots, objective=objective,
       converged=final_change(objective)<parameters$tolerance)
}

export_reference <- function(fit, output_dir) {
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  write_tsv <- function(x, name) write.table(x,file.path(output_dir,name),
      sep='\t',quote=FALSE,row.names=FALSE)
  cds <- fit$cds
  z <- t(reducedDimS(cds)); nodes <- t(reducedDimK(cds)); w <- reducedDimW(cds)
  metadata <- pData(cds)
  pt <- metadata$Pseudotime
  stopifnot(diff(range(pt))>0)
  cells <- data.frame(cell_id=rownames(metadata), metadata,
                     DDRTree1=z[,1], DDRTree2=z[,2],
                     Pseudotime01=(pt-min(pt))/diff(range(pt)), check.names=FALSE)
  projected <- t(cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_dist)
  nearest <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
  write_tsv(cells,'cells.tsv')
  write_tsv(data.frame(center_id=rownames(nodes),DDRTree1=nodes[,1],DDRTree2=nodes[,2]),'centers.tsv')
  write_tsv(igraph::as_data_frame(minSpanningTree(cds),what='edges'),'center_edges.tsv')
  write_tsv(data.frame(feature_id=rownames(cds),W1=w[,1],W2=w[,2]),'loadings.tsv')
  write_tsv(data.frame(cell_id=rownames(projected),P1=projected[,1],P2=projected[,2],
      nearest_center=rownames(nodes)[as.integer(nearest[rownames(projected),1])]),'reference_tree_projection.tsv')
  write_tsv(fit$root_candidates,'root_candidates.tsv')
  write_tsv(data.frame(iteration=seq_along(fit$objective),objective=fit$objective),'objective.tsv')
  saveRDS(cds,file.path(output_dir,'reference.rds'))
}

write_spline_basis <- function(t, grid, k, output_dir) {
  suppressPackageStartupMessages(library(mgcv))
  stopifnot(all(is.finite(t)), length(unique(t))>=k+1)
  sm <- smoothCon(s(t,bs='cr',k=k),data=data.frame(t=t),absorb.cons=TRUE)[[1]]
  B <- cbind(1,sm$X)
  G <- cbind(1,PredictMat(sm,data.frame(t=grid)))
  stopifnot(qr(B)$rank==ncol(B))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.table(B,file.path(output_dir,'design.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(G,file.path(output_dir,'prediction_design.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  invisible(list(design=B,prediction_design=G))
}

# Functions can be sourced interactively. Command-line modes use local inputs.
if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (length(args)!=4L || !args[1] %in% c('fit','basis'))
    stop('Usage: Rscript fit_reference.R fit|basis INPUT_DIR OUTPUT_DIR SETTINGS.json')
  settings <- jsonlite::fromJSON(args[4])
  input <- args[2]; output <- args[3]
  if (args[1]=='fit') {
    required <- c('ncenter','lambda_mult','sigma','gamma','max_iter','expanded_max_iter','tolerance')
    if (any(vapply(settings$trajectory[required],is.null,logical(1))) ||
        is.null(settings$baseline_stage) || is.null(settings$root_neighbors))
      stop('Replace the study-specific null settings before running')
    x <- as.matrix(read.delim(gzfile(file.path(input,'ordering.tsv.gz')),row.names=1,check.names=FALSE))
    cells <- read.delim(file.path(input,'cells.tsv'),row.names=1,check.names=FALSE)
    fit <- fit_reference(x,cells,settings$trajectory,settings$baseline_stage,settings$root_neighbors)
    export_reference(fit,output)
  } else {
    if (is.null(settings$spline_k)) stop('Set spline_k before running')
    t <- read.delim(file.path(input,'pseudotime.tsv'))$t
    grid <- read.delim(file.path(input,'prediction_grid.tsv'))$t
    write_spline_basis(t,grid,settings$spline_k,output)
  }
}
