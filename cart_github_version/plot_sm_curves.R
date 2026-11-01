# Five prespecified SMs: original-TIC response, saved native cr basis, OLS + HC3.
script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
panel_dir <- normalizePath(args[1])
source(file.path(dirname(normalizePath(script)),'plot_common.R'))
settings <- fromJSON(args[2])

dd <- file.path(panel_dir,'OriginalData')
a <- fread(file.path(dd,'abundance.tsv'));targets <- fread(file.path(dd,'target_reference.tsv'))
donors <- settings$reference_donors;colours <- setNames(c('#D58F27','#2878B5','#AC577C'),donors)
stopifnot(length(donors)==3)
curves <- list();coefficients <- list();checks <- list()
for(ds in donors) {
  bd <- file.path(dd,'basis',ds)
  axis <- fread(file.path(bd,'cell_order.tsv'));tt <- read.csv(file.path(bd,'basis_input.csv'))$t
  B <- as.matrix(read.table(file.path(bd,'basis.tsv')));G <- as.matrix(read.table(file.path(bd,'basis_grid.tsv')))
  gt <- read.csv(file.path(bd,'basis_grid.csv'))$t
  meta <- fromJSON(file.path(bd,'basis_metadata.json'))
  stopifnot(meta$basis=='cr',meta$k==settings$spline_k,!meta$penalized,max(abs(tt-axis$Pseudotime01))<1e-12)
  # Verify knots/constraints without replacing the saved design used for fitting.
  sm <- mgcv::smoothCon(mgcv::s(t,bs='cr',k=settings$spline_k),data=data.frame(t=tt),absorb.cons=TRUE)[[1]]
  error <- max(abs(sm$X-B));grid_error <- max(abs(mgcv::PredictMat(sm,data.frame(t=gt))-G))
  stopifnot(error<1e-10,grid_error<1e-10,ncol(B)==settings$spline_k-1,nrow(G)==length(gt))
  X <- cbind(1,B);P <- cbind(1,G);inv <- solve(crossprod(X));hat <- rowSums((X%*%inv)*X)
  for(target in targets$reference_row) {
    z <- a[donor==ds & target_id==target];z <- z[match(axis$cell_id,z$native_cell_id)]
    stopifnot(nrow(z)==nrow(X),all(z$native_cell_id==axis$cell_id),all(is.finite(z$abundance)),max(abs(z$Pseudotime01-tt))<1e-12)
    fit <- lm.fit(X,z$abundance);stopifnot(fit$rank==settings$spline_k)
    y <- as.numeric(P%*%fit$coefficients)
    robust <- inv%*%crossprod(X*as.numeric(fit$residuals/(1-hat)))%*%inv
    half <- 1.96*sqrt(pmax(0,rowSums((P%*%robust)*P)))
    curves[[length(curves)+1]] <- data.table(donor=ds,target_id=target,t=gt,fitted_abundance=y,HC3_halfwidth=half,lower=y-half,upper=y+half,n_cells=nrow(z))
    coefficients[[length(coefficients)+1]] <- data.table(donor=ds,target_id=target,term=c('Intercept',paste0('cr',seq_len(ncol(B)))),coefficient=fit$coefficients)
  }
  checks[[ds]] <- list(n_cells=nrow(B),knots=unname(sm$xp),rank=qr(X)$rank,basis_max_error=error,prediction_basis_max_error=grid_error,response='original TIC log1p; zeros retained')
}
d <- rbindlist(curves)
fwrite(d,file.path(dd,'fitted_curves.tsv'),sep='\t');fwrite(rbindlist(coefficients),file.path(dd,'coefficients.tsv'),sep='\t')
write_json(list(status='PASS',method='Saved cr(k=settings$spline_k) + intercept, unpenalized OLS; HC3 pointwise 95% mean-curve interval',checks=checks),file.path(dd,'spline_verification.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
plots <- lapply(seq_len(nrow(targets)),function(i) {
  z <- d[target_id==targets$reference_row[i]]
  ggplot(z,aes(t,fitted_abundance,colour=donor,fill=donor))+
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.10,colour=NA,show.legend=FALSE)+
    geom_line(linewidth=.6,show.legend=FALSE)+
    scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
    scale_x_continuous(breaks=c(0,.5,1),limits=c(0,1),expand=expansion(mult=0))+
    labs(title=targets$label[i],x=if(i==5)'Relative pseudotime' else NULL,
         y='TIC-normalized abundance')+base_theme()+theme(legend.position='none')
})
render_grid(plots,panel_dir,'figure4f',5,1,99,213,legend=setNames(colours,settings$reference_labels))
