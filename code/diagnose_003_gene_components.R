args <- commandArgs(trailingOnly=TRUE)
agg_path <- if(length(args)>=1) args[[1]] else "runs/003_extraction/sample_state_candidate_aggregates.tsv.gz"
pair_path <- if(length(args)>=2) args[[2]] else "runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv"
outdir <- if(length(args)>=3) args[[3]] else "runs/003_posthoc_component_diagnostics"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
genes <- c("AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12")
a <- read.delim(agg_path,check.names=FALSE,stringsAsFactors=FALSE)
pairs <- read.delim(pair_path,check.names=FALSE,stringsAsFactors=FALSE)
a <- a[a$cell_set=="exclude_predicted_doublets" & a$final_analysis=="Non ileal CT colonocyte",]
stopifnot(!anyDuplicated(a$sample_id))
for(g in genes) a[[g]] <- log2((a[[paste0(g,"_counts")]]+0.5)/(a$total_UMI+1)*1e6)
keep <- c("sample_id","n_cells",genes)
pre <- a[,keep]; names(pre)[-1] <- paste0("pre_",names(pre)[-1])
post <- a[,keep]; names(post)[-1] <- paste0("post_",names(post)[-1])
x <- merge(pairs,pre,by.x="pre_sample_id",by.y="sample_id",all.x=TRUE)
x <- merge(x,post,by.x="post_sample_id",by.y="sample_id",all.x=TRUE)
x$evaluable <- x$pre_n_cells>=20 & x$post_n_cells>=20
for(g in genes) x[[paste0("delta_",g)]] <- ifelse(x$evaluable,x[[paste0("post_",g)]]-x[[paste0("pre_",g)]],NA)
patient <- unique(x[,c("patient","remission","library_type","batch")])
for(g in genes) {
 one <- aggregate(x[[paste0("delta_",g)]],list(patient=x$patient),function(v) if(all(is.na(v))) NA_real_ else mean(v,na.rm=TRUE))
 names(one)[2] <- paste0("delta_",g); patient <- merge(patient,one,by="patient",all.x=TRUE)
}
exactp <- function(y,g) {
 ok <- is.finite(y) & g %in% c(0,1); y<-y[ok]; g<-g[ok]; n<-length(y); k<-sum(g==1)
 if(n<2 || k==0 || k==n) return(NA_real_)
 obs<-mean(y[g==1])-mean(y[g==0]); cmb<-combn(n,k)
 dif<-apply(cmb,2,function(ix) mean(y[ix])-mean(y[-ix]))
 mean(abs(dif)>=abs(obs)-1e-12)
}
fitone <- function(y,group,label,type) {
 ok<-is.finite(y); y<-y[ok]; group<-group[ok]; g<-as.integer(group=="Remission")
 fit<-lm(y~g); X<-model.matrix(fit); e<-resid(fit); h<-hatvalues(fit)
 br<-solve(crossprod(X)); V<-br%*%crossprod(X,X*(e^2/(1-h)^2))%*%br
 b<-coef(fit)[2]; se<-summary(fit)$coef[2,2]; seh<-sqrt(V[2,2]); df<-fit$df.residual; cr<-qt(.975,df)
 data.frame(type=type,component=label,n=length(y),remission=sum(g),nonremission=sum(1-g),
  mean_remission=mean(y[g==1]),mean_nonremission=mean(y[g==0]),estimate=b,
  lower=b-cr*se,upper=b+cr*se,p_OLS=summary(fit)$coef[2,4],
  lower_HC3=b-cr*seh,upper_HC3=b+cr*seh,p_HC3=2*pt(abs(b/seh),df,lower.tail=FALSE),
  p_exact_permutation=exactp(y,g),all_remission_above_all_nonremission=min(y[g==1])>max(y[g==0]))
}
items <- lapply(genes,function(g) fitone(patient[[paste0("delta_",g)]],patient$remission,g,"single_gene"))
for(drop in genes) {
 retained<-setdiff(genes,drop); y<-rowMeans(patient[,paste0("delta_",retained)],na.rm=FALSE)
 items[[length(items)+1]]<-fitone(y,patient$remission,paste0("drop_",drop),"leave_one_gene_out_score")
}
composite <- rowMeans(patient[,paste0("delta_",genes)],na.rm=FALSE)
items[[length(items)+1]]<-fitone(composite,patient$remission,"six_gene_mean","composite")
models <- do.call(rbind,items)
models$BH_six_genes <- NA_real_
models$BH_six_genes[models$type=="single_gene"] <- p.adjust(models$p_exact_permutation[models$type=="single_gene"],"BH")
write.table(patient,file.path(outdir,"patient_gene_deltas.tsv"),sep="	",quote=FALSE,row.names=FALSE)
write.table(models,file.path(outdir,"component_models.tsv"),sep="	",quote=FALSE,row.names=FALSE)
writeLines(c(
 "STATUS: POSTHOC_COMPONENT_DIAGNOSTIC",
 "The frozen main test is the six-gene CT pseudobulk score. Gene-level and leave-one-gene-out rows diagnose component consistency and are not additional validation tests.",
 paste0("Created UTC: ",format(Sys.time(),tz="UTC",usetz=TRUE))
),file.path(outdir,"README.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
cat("POSTHOC_COMPONENT_DIAGNOSTIC_COMPLETE")
