# Cross-language self-check of the numerical outputs; not an independent review.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
out <- 'runs/005_health_function'
pt <- fread(file.path(out,'ct_patient_scores.tsv'))
reported <- fread(file.path(out,'ct_program_models.tsv'))
checks <- list()
for(s in unique(pt$program)) {
  d <- pt[program==s]
  f <- lm(post~remission_binary+baseline,d)
  X <- model.matrix(f);bread <- solve(crossprod(X));e <- resid(f)/(1-hatvalues(f))
  V <- bread %*% crossprod(X,X*as.numeric(e^2)) %*% bread
  j <- which(colnames(X)=='remission_binary');se <- sqrt(V[j,j]);b <- coef(f)[j]
  ci <- b+c(-1,1)*qt(.975,df.residual(f))*se;pv<-2*pt(-abs(b/se),df.residual(f))
  r <- reported[program==s & model=='baseline_ANCOVA']
  er <- max(abs(c(b,se,ci,pv)-c(r$estimate,r$se,r$lower,r$upper,r$p)))
  stopifnot(er<1e-10)
  checks[[length(checks)+1]] <- data.table(check='HC3_ANCOVA',program=s,max_error=er)
}
bulk<-fread(file.path(out,'bulk_patient_scores.tsv')); controls<-fread(file.path(out,'healthy_control_units.tsv'))
hc<-fread(file.path(out,'healthy_reference_contrasts.tsv'))
for(i in which(hc$scope=='main')) {
  row<-hc[i];x<-bulk[context==row$context & program==row$program & outcome==row$outcome]
  y<-controls[series==x$series[1] & program==row$program]
  tt<-t.test(x$post,y$score)
  er<-max(abs(c(unname(tt$estimate[1]-tt$estimate[2]),tt$conf.int,tt$p.value)-c(row$estimate,row$lower,row$upper,row$p)))
  stopifnot(er<1e-10)
  checks[[length(checks)+1]]<-data.table(check='Welch_gap',program=paste(row$context,row$program,row$outcome),max_error=er)
}
primary<-c('CT_MARKERS_EXCLUDING_CANDIDATE6','GOBP_INTESTINAL_ABSORPTION','GOBP_BRUSH_BORDER_ASSEMBLY')
q<-reported[program %in% primary & model=='baseline_ANCOVA']
stopifnot(max(abs(p.adjust(q$p,'BH',n=3)-q$BH_primary_n3))<1e-12)
hm<-hc[scope=='main'];stopifnot(max(abs(p.adjust(hm$p,'BH',n=16)-hm$BH_n16))<1e-12)
membership<-fread('runs/005_preparation/program_membership.tsv')
candidate<-membership[program=='CANDIDATE6']$gene
stopifnot(nrow(membership[program!='CANDIDATE6' & gene %in% candidate])==0)
# Independently reconstruct scores and site-balanced patient summaries from counts.
counts<-as.data.table(read.delim(gzfile('runs/005_ct_extraction/sample_gene_counts.tsv.gz'),check.names=FALSE));pairs<-fread(file.path(out,'ct_eligible_site_pairs.tsv'))
for(s in unique(membership$program)) {
  genes<-membership[program==s]$gene
  y<-rowMeans(log2((as.matrix(counts[,..genes])+.5)/(counts$total_UMI+1)*1e6));names(y)<-counts$sample_id
  ss<-data.table(patient=pairs$patient,pre=y[pairs$pre_sample_id],post=y[pairs$post_sample_id])
  pp<-ss[,.(pre=mean(pre),post=mean(post)),by=patient]
  old<-pt[program==s][match(pp$patient,patient)]
  er<-max(abs(c(pp$pre-old$baseline,pp$post-old$post)))
  stopifnot(er<1e-10)
  checks[[length(checks)+1]]<-data.table(check='count_to_patient_score',program=s,max_error=er)
}
qc<-fread(file.path(out,'ct_program_detection.tsv'))
summ<-qc[,.(members=.N,never_detected=sum(n_eligible_samples_detected==0),median_detected_cell_fraction=median(fraction_CT_cells_detected),members_seen_in_all_66=sum(n_eligible_samples_detected==66)),by=program]
fwrite(summ,file.path(out,'program_detection_summary.tsv'),sep='\t')
fwrite(rbindlist(checks),file.path(out,'cross_language_checks.tsv'),sep='\t')
write_json(list(status='PASS',review_type='cross-language self-check',n_checks=length(checks),max_error=max(rbindlist(checks)$max_error),BHadjustments='n3 and n16 verified',candidate_excluded_from_references=TRUE,R=R.version.string),file.path(out,'verification.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
print(summ)
