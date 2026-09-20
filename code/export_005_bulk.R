options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root <- normalizePath('.')
out <- file.path(root,'runs/005_bulk_export')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'COMPLETE.json'))) stop('005 bulk export already complete')
plan <- read_json('planning/analysis_plan_005_health_function.json',simplifyVector=TRUE)
stopifnot(plan$status=='FROZEN')
d <- readRDS('runs/001D_annotation/discovery_preprocessed_canonical.rds')
b <- readRDS('runs/002B_GSE73661_support/support_scores.rds')
reg <- fread('runs/001A_metadata/sample_registry_v1.tsv')
sets <- plan$programs
scores <- list(); coverage <- list()
for(series in c('GSE92415','GSE73661')) {
  g <- if(series=='GSE92415') d$gene_expression else b$gene_expression
  for(prog in names(sets)) {
    hit <- intersect(sets[[prog]],rownames(g))
    coverage[[length(coverage)+1]] <- data.table(series=series,program=prog,n_members=length(sets[[prog]]),n_measured=length(hit),coverage=length(hit)/length(sets[[prog]]),missing=paste(setdiff(sets[[prog]],hit),collapse=';'))
    if(prog %in% c('CANDIDATE6','COMMON194_INFLAMMATION')) stopifnot(length(hit)==length(sets[[prog]]))
    if(length(hit)>=5 && length(hit)/length(sets[[prog]])>=.8) {
      scores[[length(scores)+1]] <- data.table(series=series,gsm=colnames(g),program=prog,score=colMeans(g[hit,,drop=FALSE]))
    }
  }
}
s <- rbindlist(scores)
fwrite(s,file.path(out,'sample_scores.tsv'),sep='\t')
fwrite(rbindlist(coverage),file.path(out,'program_coverage.tsv'),sep='\t')
fwrite(reg[series %in% c('GSE92415','GSE73661')],file.path(out,'sample_registry.tsv'),sep='\t')
check <- max(abs(s[series=='GSE73661' & program=='CANDIDATE6']$score-b$sample_scores['MATURE_COLONOCYTE6_EXPLORATORY',s[series=='GSE73661' & program=='CANDIDATE6']$gsm]))
stopifnot(check<1e-10)
write_json(list(status='COMPLETE',score_max_error_vs_002B=check,R=R.version.string,utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))

