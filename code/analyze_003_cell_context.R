args <- commandArgs(trailingOnly=TRUE)
agg_path <- if(length(args)>=1) args[[1]] else "runs/003_extraction/sample_state_candidate_aggregates.tsv.gz"
pair_path <- if(length(args)>=2) args[[2]] else "runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv"
outdir <- if(length(args)>=3) args[[3]] else "runs/003_cell_context"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
genes <- c("AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12")
ag <- read.delim(agg_path, check.names=FALSE, stringsAsFactors=FALSE)
pairs <- read.delim(pair_path, check.names=FALSE, stringsAsFactors=FALSE)
required <- c("cell_set","sample_id","Patient","Site","Treatment","Remission_status",
              "LibraryType","Batch","final_analysis","major","n_cells","total_UMI",
              paste0(genes,"_counts"),paste0(genes,"_norm10k_sum"))
stopifnot(all(required %in% names(ag)))

sum_sample <- function(z, prefix) {
  if(nrow(z)==0) return(NULL)
  d <- data.frame(sample_id=unique(z$sample_id))
  stopifnot(nrow(d)==1)
  d[[paste0(prefix,"_n")]] <- sum(z$n_cells)
  d[[paste0(prefix,"_UMI")]] <- sum(z$total_UMI)
  for(g in genes) {
    d[[paste0(prefix,"_",g,"_counts")]] <- sum(z[[paste0(g,"_counts")]])
    d[[paste0(prefix,"_",g,"_normsum")]] <- sum(z[[paste0(g,"_norm10k_sum")]])
  }
  d
}
make_sample <- function(a) {
  meta <- unique(a[,c("sample_id","Patient","Site","Treatment","Remission_status","LibraryType","Batch")])
  stopifnot(!anyDuplicated(meta$sample_id))
  allx <- do.call(rbind,lapply(split(a,a$sample_id),sum_sample,prefix="all"))
  epi <- a[a$major=="Non_ileal_epithelium",]
  epix <- do.call(rbind,lapply(split(epi,epi$sample_id),sum_sample,prefix="epi"))
  ct <- a[a$final_analysis=="Non ileal CT colonocyte",]
  ctx <- do.call(rbind,lapply(split(ct,ct$sample_id),sum_sample,prefix="ct"))
  out <- merge(meta,allx,by="sample_id",all.x=TRUE)
  out <- merge(out,epix,by="sample_id",all.x=TRUE)
  out <- merge(out,ctx,by="sample_id",all.x=TRUE)
  countcols <- grep("_(n|UMI|counts|normsum)$",names(out),value=TRUE)
  out[countcols][is.na(out[countcols])] <- 0
  out$ct_logit_epi <- ifelse(out$epi_n>0,log((out$ct_n+0.5)/(out$epi_n-out$ct_n+0.5)),NA)
  out$ct_logit_all <- ifelse(out$all_n>0,log((out$ct_n+0.5)/(out$all_n-out$ct_n+0.5)),NA)
  ct_scores <- sapply(genes,function(g) log2((out[[paste0("ct_",g,"_counts")]]+0.5)/(out$ct_UMI+1)*1e6))
  epi_scores <- sapply(genes,function(g) log2((out[[paste0("epi_",g,"_counts")]]+0.5)/(out$epi_UMI+1)*1e6))
  out$ct_state_score <- rowMeans(ct_scores)
  out$epi_score <- rowMeans(epi_scores)
  out$ct_linear <- rowMeans(sapply(genes,function(g) out[[paste0("ct_",g,"_normsum")]]/pmax(out$ct_n,1)))
  other_n <- out$epi_n-out$ct_n
  out$other_n <- other_n
  out$other_linear <- rowMeans(sapply(genes,function(g)
    (out[[paste0("epi_",g,"_normsum")]]-out[[paste0("ct_",g,"_normsum")]])/pmax(other_n,1)))
  out$epi_linear <- rowMeans(sapply(genes,function(g) out[[paste0("epi_",g,"_normsum")]]/pmax(out$epi_n,1)))
  out
}
pair_samples <- function(s, pairtab, threshold) {
  pre <- s; names(pre)[names(pre)!="sample_id"] <- paste0("pre_",names(pre)[names(pre)!="sample_id"])
  post <- s; names(post)[names(post)!="sample_id"] <- paste0("post_",names(post)[names(post)!="sample_id"])
  x <- merge(pairtab,pre,by.x="pre_sample_id",by.y="sample_id",all.x=TRUE)
  x <- merge(x,post,by.x="post_sample_id",by.y="sample_id",all.x=TRUE)
  x$threshold <- threshold
  x$delta_ct_logit_epi <- x$post_ct_logit_epi-x$pre_ct_logit_epi
  x$delta_ct_logit_all <- x$post_ct_logit_all-x$pre_ct_logit_all
  eligible <- x$pre_ct_n>=threshold & x$post_ct_n>=threshold
  x$state_evaluable <- eligible
  x$delta_ct_state_score <- ifelse(eligible,x$post_ct_state_score-x$pre_ct_state_score,NA)
  x$delta_epi_score <- x$post_epi_score-x$pre_epi_score
  decomp_ok <- eligible & x$pre_other_n>0 & x$post_other_n>0
  p0 <- x$pre_ct_n/x$pre_epi_n; p1 <- x$post_ct_n/x$post_epi_n
  x$composition_contribution <- ifelse(decomp_ok,(p1-p0)*((x$post_ct_linear+x$pre_ct_linear)-
      (x$post_other_linear+x$pre_other_linear))/2,NA)
  x$within_state_contribution <- ifelse(decomp_ok,((p1+p0)/2)*(x$post_ct_linear-x$pre_ct_linear)+
      ((2-p1-p0)/2)*(x$post_other_linear-x$pre_other_linear),NA)
  x$total_epi_linear_change <- ifelse(decomp_ok,x$post_epi_linear-x$pre_epi_linear,NA)
  x$decomposition_error <- x$total_epi_linear_change-x$composition_contribution-x$within_state_contribution
  x
}
patient_means <- function(x) {
  metrics <- c("delta_ct_logit_epi","delta_ct_logit_all","delta_ct_state_score","delta_epi_score",
               "composition_contribution","within_state_contribution","total_epi_linear_change",
               "decomposition_error")
  z <- aggregate(x[,metrics],by=list(patient=x$patient,remission=x$remission,
                 library_type=x$library_type,batch=x$batch),FUN=function(v) if(all(is.na(v))) NA else mean(v,na.rm=TRUE))
  nsite <- aggregate(rep(1,nrow(x)),by=list(patient=x$patient),FUN=sum)
  names(nsite)[2] <- "n_site_pairs"
  merge(z,nsite,by="patient")
}
exact_label_p <- function(y,g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]; g <- g[ok]
  n <- length(y); k <- sum(g==1)
  if(n<2 || k==0 || k==n) return(NA_real_)
  obs <- mean(y[g==1])-mean(y[g==0])
  cmb <- combn(n,k)
  dif <- apply(cmb,2,function(ix) mean(y[ix])-mean(y[-ix]))
  mean(abs(dif)>=abs(obs)-1e-12)
}
fit_group <- function(d, metric, branch) {
  ok <- is.finite(d[[metric]]) & d$remission %in% c("Remission","Non_Remission")
  z <- d[ok,]; z$g <- as.integer(z$remission=="Remission")
  if(nrow(z)<4 || length(unique(z$g))<2) return(data.frame(branch=branch,metric=metric,n=nrow(z),
      remission=sum(z$g),nonremission=sum(1-z$g),estimate=NA,se=NA,lower=NA,upper=NA,p_OLS=NA,
      se_HC3=NA,lower_HC3=NA,upper_HC3=NA,p_HC3=NA,p_exact_permutation=NA))
  fit <- lm(z[[metric]]~z$g); X <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  bread <- solve(crossprod(X)); meat <- crossprod(X,X*(e^2/(1-h)^2)); V <- bread%*%meat%*%bread
  b <- coef(fit)[2]; se <- summary(fit)$coef[2,2]; seh <- sqrt(V[2,2]); df <- fit$df.residual
  crit <- qt(.975,df)
  data.frame(branch=branch,metric=metric,n=nrow(z),remission=sum(z$g),nonremission=sum(1-z$g),
    mean_remission=mean(z[[metric]][z$g==1]),mean_nonremission=mean(z[[metric]][z$g==0]),
    estimate=b,se=se,lower=b-crit*se,upper=b+crit*se,p_OLS=summary(fit)$coef[2,4],
    se_HC3=seh,lower_HC3=b-crit*seh,upper_HC3=b+crit*seh,
    p_HC3=2*pt(abs(b/seh),df=df,lower.tail=FALSE),
    p_exact_permutation=exact_label_p(z[[metric]],z$g))
}
sign_p <- function(y) {
 y <- y[is.finite(y)]; n <- length(y)
 if(n==0 || n>22) return(NA_real_)
 signs <- as.matrix(expand.grid(rep(list(c(-1,1)),n)))
 obs <- abs(mean(y)); mean(abs(as.vector(signs%*%y/n))>=obs-1e-12)
}
fit_change <- function(d,metric,branch) {
 y=d[[metric]]; y=y[is.finite(y)]; n=length(y)
 if(n<2) return(data.frame(branch=branch,metric=metric,n=n,mean_change=NA,se=NA,lower=NA,upper=NA,p_t=NA,p_sign_permutation=NA))
 se=sd(y)/sqrt(n); crit=qt(.975,n-1); est=mean(y)
 data.frame(branch=branch,metric=metric,n=n,mean_change=est,se=se,lower=est-crit*se,upper=est+crit*se,
            p_t=2*pt(abs(est/se),df=n-1,lower.tail=FALSE),p_sign_permutation=if(branch=="exclude_predicted_doublets::all_fixed::threshold20") sign_p(y) else NA_real_)
}
all_sites <- list(); all_patients <- list(); models <- list(); changes <- list()
for(cellset in c("exclude_predicted_doublets","author_annotated")) {
  a <- ag[ag$cell_set==cellset,]; s <- make_sample(a)
  for(thr in c(10,20,50)) {
    for(scope in c("all_fixed","v31_only","author_both")) {
      pr <- pairs
      if(scope=="v31_only") pr <- pr[pr$library_type=="3 prime v3.1",]
      if(scope=="author_both") pr <- pr[tolower(as.character(pr$both_in_author_paired_list))=="true",]
      site <- pair_samples(s,pr,thr); site$cell_set <- cellset; site$scope <- scope
      patient <- patient_means(site)
      branch <- paste(cellset,scope,paste0("threshold",thr),sep="::")
      cat("BRANCH",branch,"site",nrow(site),"patient",nrow(patient),"\n")
      if(nrow(patient)==0) stop("No patient rows in ",branch)
      patient$cell_set <- cellset; patient$scope <- scope; patient$threshold <- thr
      for(metric in c("delta_ct_logit_epi","delta_ct_logit_all","delta_ct_state_score","delta_epi_score",
                      "composition_contribution","within_state_contribution","total_epi_linear_change")) {
        models[[length(models)+1]] <- fit_group(patient,metric,branch)
        changes[[length(changes)+1]] <- fit_change(patient,metric,branch)
      }
      all_sites[[length(all_sites)+1]] <- site
      all_patients[[length(all_patients)+1]] <- patient
    }
  }
}
site_all <- do.call(rbind,all_sites); patient_all <- do.call(rbind,all_patients)
model_all <- do.call(rbind,models); change_all <- do.call(rbind,changes)
primary <- model_all$branch=="exclude_predicted_doublets::all_fixed::threshold20" &
           model_all$metric %in% c("delta_ct_logit_epi","delta_ct_state_score")
model_all$BH_two_main <- NA_real_
model_all$BH_two_main[primary] <- p.adjust(model_all$p_exact_permutation[primary],method="BH")
write.table(site_all,file.path(outdir,"site_pair_metrics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(patient_all,file.path(outdir,"patient_metrics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(model_all,file.path(outdir,"patient_group_models.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(change_all,file.path(outdir,"all_patient_change_models.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# Fixed-family non-ileal epithelial state mapping, main cell set and threshold20.
a <- ag[ag$cell_set=="exclude_predicted_doublets" & ag$major=="Non_ileal_epithelium",]
for(g in genes) a[[paste0(g,"_logCPM")]] <- log2((a[[paste0(g,"_counts")]]+0.5)/(a$total_UMI+1)*1e6)
a$state_score <- rowMeans(a[,paste0(genes,"_logCPM")])
pre <- a[,c("sample_id","final_analysis","n_cells","state_score")]
names(pre)[3:4] <- c("pre_n","pre_score")
post <- pre; names(post)[3:4] <- c("post_n","post_score")
st <- merge(pairs,pre,by.x="pre_sample_id",by.y="sample_id",all=FALSE)
st <- merge(st,post,by.x=c("post_sample_id","final_analysis"),by.y=c("sample_id","final_analysis"),all=FALSE)
st <- st[st$pre_n>=20 & st$post_n>=20,]
st$delta_state_score <- st$post_score-st$pre_score
sp <- aggregate(delta_state_score~patient+remission+final_analysis,data=st,FUN=mean)
states <- sort(unique(a$final_analysis))
state_models <- do.call(rbind,lapply(states,function(state) fit_group(sp[sp$final_analysis==state,],
                                                                     "delta_state_score",state)))
state_models$BH_fixed_state_family <- p.adjust(state_models$p_exact_permutation,method="BH")
write.table(state_models,file.path(outdir,"fixed_epithelial_state_models.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

maxerr <- max(abs(site_all$decomposition_error),na.rm=TRUE)
summary <- list(status="COMPLETE",role="Exploratory cell-context discrimination; no confirmatory primary tests",
 n_fixed_patients=length(unique(pairs$patient)),n_remission=length(unique(pairs$patient[pairs$remission=="Remission"])),
 primary_branch="exclude_predicted_doublets::all_fixed::threshold20",
 primary_metrics=c("delta_ct_logit_epi","delta_ct_state_score"),
 max_absolute_decomposition_error=maxerr,
 chemistry_boundary="All remitters are v3.1; v3.1-only branch is mandatory sensitivity, not a complete solution",
 cross_cohort_boundary="Adalimumab clinical-remission context cannot settle GSE73661 IFX endoscopic construct overlap")
jsonlite::write_json(summary,file.path(outdir,"analysis_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
cat("STAGE003_CELL_CONTEXT_COMPLETE\n")