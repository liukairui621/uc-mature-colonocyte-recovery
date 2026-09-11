args <- commandArgs(trailingOnly=TRUE)
if(length(args)<6) stop("args: sample_id npz metadata var outdir seed")
sample_id<-args[1]; npz<-args[2]; meta_path<-args[3]; var_path<-args[4]; outdir<-args[5]; seed<-as.integer(args[6])
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
status_path<-file.path(outdir,paste0(sample_id,".status.json"))
started<-format(Sys.time(),tz="UTC",usetz=TRUE)
status<-list(status="RUNNING",sample_id=sample_id,started_utc=started,seed=seed)
jsonlite::write_json(status,status_path,auto_unbox=TRUE,pretty=TRUE)
tryCatch({
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(reticulate))
  use_python("/usr/bin/python3",required=TRUE)
  sp<-import("scipy.sparse",convert=TRUE)
  cell_by_gene<-sp$load_npz(npz)
  counts<-as(t(cell_by_gene),"dgCMatrix")
  md<-read.delim(meta_path,check.names=FALSE,stringsAsFactors=FALSE)
  var<-read.delim(var_path,check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(ncol(counts)==nrow(md),nrow(counts)==nrow(var))
  genes<-c("AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12")
  gi<-match(genes,var$gene_symbol)
  if(anyNA(gi)) stop("Candidate gene missing")
  z<-md$final_analysis
  if(anyNA(z)||length(unique(z))<2) stop("Invalid final-analysis clusters")
  log_path<-file.path(outdir,paste0(sample_id,".decontx.log"))
  res<-celda::decontX(counts,z=z,maxIter=500,delta=c(10,10),estimateDelta=TRUE,
                      convergence=0.001,iterLogLik=10,varGenes=5000,
                      seed=seed,logfile=log_path,verbose=TRUE)
  dc<-res$decontXcounts
  contamination<-as.numeric(res$contamination)
  lib<-Matrix::colSums(dc)
  states<-sort(unique(z))
  agg<-do.call(rbind,lapply(states,function(st){
    ix<-which(z==st)
    one<-data.frame(sample_id=sample_id,Patient=md$Patient[1],Site=md$Site[1],
      Treatment=md$Treatment[1],Remission_status=md$Remission_status[1],
      LibraryType=md$LibraryType[1],Batch=md$Batch[1],final_analysis=st,
      major=md$major[ix[1]],n_cells=length(ix),decontx_total_UMI=sum(lib[ix]),
      contamination_mean=mean(contamination[ix]),contamination_median=median(contamination[ix]),
      contamination_q90=unname(quantile(contamination[ix],.9)),stringsAsFactors=FALSE)
    for(j in seq_along(genes)) one[[paste0(genes[j],"_decontx_counts")]]<-sum(dc[gi[j],ix])
    one
  }))
  write.table(agg,file.path(outdir,paste0(sample_id,".sample_state_decontx.tsv")),
              sep="	",quote=FALSE,row.names=FALSE)
  status<-list(status="COMPLETE",sample_id=sample_id,started_utc=started,
    completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),seed=seed,
    n_cells=ncol(counts),n_genes=nrow(counts),n_clusters=length(unique(z)),
    celda_version=as.character(packageVersion("celda")),
    matrix_version=as.character(packageVersion("Matrix")),
    contamination_mean=mean(contamination),contamination_median=median(contamination),
    contamination_q90=unname(quantile(contamination,.9)),
    corrected_total_UMI=sum(lib),raw_total_UMI=sum(counts))
  jsonlite::write_json(status,status_path,auto_unbox=TRUE,pretty=TRUE,digits=16)
},error=function(e){
  status<-list(status="FAILED",sample_id=sample_id,started_utc=started,
    failed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),seed=seed,error=conditionMessage(e))
  jsonlite::write_json(status,status_path,auto_unbox=TRUE,pretty=TRUE)
  stop(e)
})
cat("DECONTX_SAMPLE_COMPLETE",sample_id)
