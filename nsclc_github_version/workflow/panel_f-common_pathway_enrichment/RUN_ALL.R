a<-grep('^--file=',commandArgs(),value=TRUE);base<-dirname(normalizePath(sub('^--file=','',a[1])))
status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(base,'scripts','01_run_KEGG_ORA.R')))
if(status!=0)stop('Panel generation failed')
