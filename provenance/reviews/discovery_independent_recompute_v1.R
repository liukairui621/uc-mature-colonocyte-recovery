options(stringsAsFactors=FALSE, width=150)
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(limma)})
root <- '/root/projects/UC_Treatment_Recovery'
out <- file.path(root,'provenance/reviews')
obj <- readRDS(file.path(root,'runs/001A_qc/discovery_preprocessed_v1.rds'))
# Independently parse matrix body with data.table; no call to production readers.
lines <- readLines(gzfile(file.path(root,'inputs/discovery/GSE92415_series_matrix.relay.txt.gz')),warn=FALSE)
begin <- which(lines=='!series_matrix_table_begin');end <- which(lines=='!series_matrix_table_end')
rawtab <- fread(text=paste(lines[(begin+1):(end-1)],collapse='\n'),check.names=FALSE)
raw <- as.matrix(rawtab[,-1]);rownames(raw)<-rawtab[[1]]
source_rma <- any(grepl('normalized by RMA using ArrayStudio version 9', lines, fixed=TRUE))
stopifnot(identical(dim(raw),dim(obj$probe_expression)),identical(dimnames(raw),dimnames(obj$probe_expression)))
raw_error <- max(abs(raw-obj$probe_expression))
rm(rawtab,lines);gc()
# Parse downloaded GPL table and collapse all retained probes by row sums/counts,
# independently of limma::avereps used by production.
alines <- readLines(gzfile(file.path(root,'inputs/annotation/GPL13158.annot.gz')),warn=FALSE)
h <- grep('^ID\t',alines)[1]
alines <- alines[h:length(alines)];alines <- alines[!startsWith(alines,'!')]
a <- fread(text=paste(alines,collapse='\n'),quote='',fill=TRUE,sep='\t',check.names=FALSE)
symbol <- a[['Gene symbol']][match(rownames(raw),a[['ID']])]
good <- !is.na(symbol) & nchar(symbol)>0 & !grepl('///|//|;|\\|',symbol) & !startsWith(rownames(raw),'AFFX')
counts <- table(symbol[good])
summed <- rowsum(raw[good,,drop=FALSE],symbol[good],reorder=TRUE)
genes <- summed / as.numeric(counts[rownames(summed)])
genes <- genes[rownames(obj$gene_expression),,drop=FALSE]
collapse_error <- max(abs(genes-obj$gene_expression))
map <- fread(file.path(root,'runs/001A_qc/probe_gene_mapping.tsv'))
map_identity <- identical(as.character(map$probe),rownames(raw)) && identical(as.logical(map$retained),good)
placeholder_symbols <- intersect(rownames(genes),c('---','NA','N/A','null','NULL','-'))
# Compare model pairs against independently parsed metadata audit, not just pipeline TSV.
meta_review <- fromJSON(file.path(out,'metadata_independent_v1.json'),simplifyVector=FALSE)
expected <- setNames(meta_review$exact_pairs,vapply(meta_review$exact_pairs,`[[`,character(1),'subject'))
pairs <- fread(file.path(root,'runs/001A_metadata/discovery_pairs_v1.tsv'))
identity_checks <- vapply(seq_len(nrow(pairs)),function(i){
 p<-pairs[i];r<-expected[[p$subject]]
 identical(p$baseline_gsm,r$baseline$gsm) && identical(p$post_gsm,r$week6$gsm) && identical(p$arm,r$baseline$treatment) && identical(p$response,r$baseline$wk6response)
},logical(1))
delta <- genes[,pairs$post_gsm,drop=FALSE]-genes[,pairs$baseline_gsm,drop=FALSE];colnames(delta)<-pairs$subject
pair_error <- max(abs(delta-obj$delta))
qc_summary <- fromJSON(file.path(root,'runs/001A_qc/qc_summary.json'))
qc_counts_match <- nrow(raw)==qc_summary$probes && ncol(raw)==qc_summary$arrays && sum(good)==qc_summary$mapped_single_gene_probes && nrow(genes)==qc_summary$genes && nrow(pairs)==qc_summary$pairs
# Deterministic genes spanning alphabetically ordered universe, chosen without P values.
sel <- unique(round(seq(1,nrow(genes),length.out=12)))
selected <- rownames(genes)[sel]
checks <- list(); recomputed <- list(); gene_count_summaries<-list()
manual_bh <- function(p){o<-order(p);adj<-pmin(1,rev(cummin(rev(p[o]*length(p)/seq_along(p)))));z<-numeric(length(p));z[o]<-adj;z}
for(branch in c('primary_paired','barcode_sensitivity')){
 bdir <- file.path(root,'runs/001B_discovery',branch)
 d <- fread(file.path(bdir,'main_design.tsv'));z <- fread(file.path(bdir,'response_design.tsv'))
 X<-as.matrix(d[,-1]);Z<-as.matrix(z[,-1]);stopifnot(identical(d$subject,pairs$subject),identical(z$subject,pairs$subject))
 expectedX <- cbind(Placebo=as.integer(pairs$arm=='Placebo'),GLM=as.integer(pairs$arm=='golimumab'))
 expectedZ <- sapply(c('Placebo_No','Placebo_Yes','golimumab_No','golimumab_Yes'),function(k)as.integer(paste(pairs$arm,pairs$response,sep='_')==k))
 if(branch=='barcode_sensitivity'){
  bm<-fread(file.path(out,'metadata_independent_v1_gsm_barcode.tsv'))
  bl<-bm$inferred_cel_barcode[match(pairs$baseline_gsm,bm$gsm)];po<-bm$inferred_cel_barcode[match(pairs$post_gsm,bm$gsm)];lev<-sort(unique(c(bl,po)))
  B<-sapply(lev[-1],function(k)as.integer(po==k)-as.integer(bl==k));expectedX<-cbind(expectedX,B);expectedZ<-cbind(expectedZ,B)
 }
 stopifnot(max(abs(X-expectedX))==0,max(abs(Z-expectedZ))==0)
 fits<-readRDS(file.path(bdir,'limma_fits.rds'))
 kinds <- list(main=list(design=X,fit=fits$main),response=list(design=Z,fit=fits$response))
 for(kind in names(kinds)){
  M<-kinds[[kind]]$design;fit<-kinds[[kind]]$fit
  beta<-qr.solve(M,t(delta[selected,,drop=FALSE]));resid<-t(delta[selected,,drop=FALSE])-M%*%beta
  sigma<-sqrt(colSums(resid^2)/(nrow(M)-ncol(M)))
  V<-solve(crossprod(M))
  for(co in colnames(fit$coefficients)){
   cvec<-setNames(rep(0,ncol(M)),colnames(M))
   if(co=='GLM_change')cvec['GLM']<-1
   if(co=='Placebo_change')cvec['Placebo']<-1
   if(co=='GLM_minus_Placebo')cvec[c('GLM','Placebo')]<-c(1,-1)
   if(co=='GLM_response_difference')cvec[c('G_Yes','G_No')]<-c(1,-1)
   if(co=='Placebo_response_difference')cvec[c('P_Yes','P_No')]<-c(1,-1)
   if(co=='response_difference_interaction')cvec[c('G_Yes','G_No','P_Yes','P_No')]<-c(1,-1,-1,1)
   estimate<-as.numeric(crossprod(cvec,beta));unscaled<-sqrt(as.numeric(crossprod(cvec,V%*%cvec)))
   ix<-match(selected,rownames(fit$coefficients));tt<-fread(file.path(bdir,paste0(co,'.tsv')));j<-match(selected,tt$gene)
   moderated_se<-unscaled*sqrt(fit$s2.post[ix]);tval<-estimate/moderated_se
   pval<-2*pt(-abs(tval),df=fit$df.total[ix]);ci<-qt(.975,df=fit$df.total[ix])*moderated_se
   results<-data.frame(branch=branch,contrast=co,gene=selected,estimate_independent=estimate,estimate_exported=tt$logFC[j],residual_sd_independent=sigma,residual_sd_fit=fit$sigma[ix],p_from_stored_moderated_variance=pval,p_exported=tt$P.Value[j])
   recomputed[[paste(branch,co)]]<-results
   checks[[paste(branch,co)]]<-list(design_rank=qr(M)$rank,design_columns=ncol(M),residual_df=nrow(M)-ncol(M),coefficient_max_abs_error=max(abs(estimate-tt$logFC[j])),sigma_max_abs_error=max(abs(sigma-fit$sigma[ix])),unscaled_se_max_abs_error=max(abs(unscaled-fit$stdev.unscaled[ix,co])),moderated_t_max_abs_error=max(abs(tval-tt$t[j])),moderated_p_max_abs_error=max(abs(pval-tt$P.Value[j])),ci_max_abs_error=max(abs(c(estimate-ci-tt$CI.L[j],estimate+ci-tt$CI.R[j]))),BH_max_abs_error=max(abs(manual_bh(tt$P.Value)-tt$adj.P.Val)),AveExpr_is_overall_mean_delta=max(abs(tt$AveExpr-rowMeans(delta)[match(tt$gene,rownames(delta))]))<1e-12)
   gene_count_summaries[[paste(branch,co)]]<-list(tested_genes=nrow(tt),fdr05=sum(tt$adj.P.Val<.05),fdr05_up=sum(tt$adj.P.Val<.05&tt$logFC>0),fdr05_down=sum(tt$adj.P.Val<.05&tt$logFC<0))
  }
 }
}
# Audit all exported CAMERA FDR, counts and coefficient summaries. Re-run one full
# Hallmark contrast under identical parameters, not alternate outcome-seeking models.
sets<-readRDS(file.path(root,'inputs/annotation/gene_sets_001B.rds'))
path_checks<-list();path_counts<-list()
for(branch in c('primary_paired','barcode_sensitivity')){
 bdir<-file.path(root,'runs/001B_discovery',branch)
 for(lib in c('Hallmark','GO_BP')){
  alltab<-fread(file.path(bdir,paste0('camera_',lib,'_all.tsv')))
  global_error<-max(abs(manual_bh(alltab$PValue)-alltab$FDR_across_four_contrasts))
  for(co in unique(alltab$contrast)){
   t<-alltab[contrast==co];expr<-fread(file.path(bdir,paste0(co,'.tsv')))
   ns<-vapply(t$set,function(k)sum(rownames(delta)%in%sets[[lib]][[k]]),integer(1))
   means<-vapply(t$set,function(k)mean(expr$logFC[expr$gene%in%sets[[lib]][[k]]]),numeric(1))
   path_checks[[paste(branch,lib,co)]]<-list(sets=nrow(t),BH_max_abs_error=max(abs(manual_bh(t$PValue)-t$FDR)),four_contrast_BH_max_abs_error=global_error,set_sizes_match=identical(as.integer(ns),as.integer(t$NGenes)),mean_coefficient_max_abs_error=max(abs(means-t$mean_log2_change)),fdr05=sum(t$FDR<.05),four_contrast_fdr05=sum(t$FDR_across_four_contrasts<.05))
   path_counts[[paste(branch,lib,co)]]<-sum(t$FDR<.05)
  }
 }
}
ix<-lapply(sets$Hallmark,function(s)which(rownames(delta)%in%s));ix<-ix[lengths(ix)>=15&lengths(ix)<=500]
d<-fread(file.path(root,'runs/001B_discovery/primary_paired/main_design.tsv'));X<-as.matrix(d[,-1]);co<-setNames(c(0,1),colnames(X))
cam<-camera(delta,index=ix,design=X,contrast=co,inter.gene.cor=.01,sort=TRUE,use.ranks=FALSE)
pub<-fread(file.path(root,'runs/001B_discovery/primary_paired/camera_Hallmark_GLM_change.tsv'));j<-match(rownames(cam),pub$set)
camcheck<-list(test='Full Hallmark primary GLM_change repeated at identical parameters',P_max_abs_error=max(abs(cam$PValue-pub$PValue[j])),FDR_max_abs_error=max(abs(cam$FDR-pub$FDR[j])),direction_identical=identical(as.character(cam$Direction),pub$Direction[j]))
# Recompute counts recorded in summaries, without assuming JSON values are correct.
ms<-fromJSON(file.path(root,'runs/001B_discovery/discovery_model_summary.json'),simplifyVector=FALSE);count_mismatches<-list()
for(branch in c('primary_paired','barcode_sensitivity')){
 original<-if(branch=='primary_paired')ms$primary$results else ms$sensitivity$results
 for(co in names(original))for(field in c('tested_genes','fdr05','fdr05_up','fdr05_down')){
  if(original[[co]][[field]]!=gene_count_summaries[[paste(branch,co)]][[field]])count_mismatches[[paste(branch,co,field)]]<-TRUE
 }
}
fwrite(rbindlist(recomputed),file.path(out,'discovery_independent_coefficients_v1.tsv'),sep='\t')
report<-list(status='COMPLETE',scope='Independent discovery-only numerical audit; no external validation expression read',raw_matrix=list(probes=nrow(raw),samples=ncol(raw),RMA_statement_in_source=source_rma,maximum_raw_vs_RDS_difference=raw_error,range=range(raw),any_nonfinite=any(!is.finite(raw))),mapping=list(retained_probes=sum(good),genes=nrow(genes),retention_identical=map_identity,collapse_max_abs_error=collapse_error,placeholder_symbols=placeholder_symbols),pairing=list(subjects=nrow(pairs),all_pairs_match_independent_header_audit=all(identity_checks),delta_max_abs_error=pair_error),qc_count_summary_identical=qc_counts_match,selected_genes=selected,coefficient_checks=checks,gene_count_summaries=gene_count_summaries,gene_count_summary_mismatches=count_mismatches,pathway_checks=path_checks,CAMERA_repeat=camcheck,session=capture.output(sessionInfo()))
write_json(report,file.path(out,'discovery_independent_numeric_v1.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat(toJSON(list(raw_max_error=raw_error,collapse_max_error=collapse_error,delta_max_error=pair_error,qc_count_summary_identical=qc_counts_match,coefficient_max_error=max(vapply(checks,`[[`,numeric(1),'coefficient_max_abs_error')),BH_max_error=max(vapply(checks,`[[`,numeric(1),'BH_max_abs_error')),summary_mismatch_count=length(count_mismatches),CAMERA_repeat=camcheck),auto_unbox=TRUE,pretty=TRUE),'\n')
