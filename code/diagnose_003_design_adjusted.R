outdir <- "runs/003_posthoc_design_diagnostics"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
p <- read.delim("runs/003_cell_context/patient_metrics.tsv",check.names=FALSE,stringsAsFactors=FALSE)
z <- p[p$cell_set=="exclude_predicted_doublets" & p$scope=="v31_only" & p$threshold==20,]
fit_hc3 <- function(metric) {
 d <- z[is.finite(z[[metric]]),]
 d$g <- as.integer(d$remission=="Remission"); d$batch <- factor(d$batch)
 fit <- lm(d[[metric]]~g+batch,data=d); X<-model.matrix(fit); e<-resid(fit); h<-hatvalues(fit)
 bread<-solve(crossprod(X)); V<-bread%*%crossprod(X,X*(e^2/(1-h)^2))%*%bread
 b<-coef(fit)["g"]; se<-summary(fit)$coef["g","Std. Error"]; seh<-sqrt(V["g","g"]); df<-fit$df.residual; cr<-qt(.975,df)
 data.frame(metric=metric,n=nrow(d),remission=sum(d$g),nonremission=sum(1-d$g),
            parameters=length(coef(fit)),residual_df=df,estimate=b,
            lower=b-cr*se,upper=b+cr*se,p_OLS=summary(fit)$coef["g","Pr(>|t|)"],
            lower_HC3=b-cr*seh,upper_HC3=b+cr*seh,p_HC3=2*pt(abs(b/seh),df,lower.tail=FALSE))
}
res <- do.call(rbind,lapply(c("delta_ct_logit_epi","delta_ct_state_score"),fit_hc3))
write.table(res,file.path(outdir,"v31_batch_adjusted_models.tsv"),sep="	",quote=FALSE,row.names=FALSE)
write.table(as.data.frame.matrix(table(z$remission,z$batch)),file.path(outdir,"v31_patient_batch_crosswalk.tsv"),sep="	",quote=FALSE)
writeLines(c("STATUS: POSTHOC_DESCRIPTIVE","Model: outcome change ~ remission + released Batch, restricted to v3.1.",
 "This cannot remove unmeasured time or batch-associated biology and is not a new validation test."),
 file.path(outdir,"README.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
cat("DESIGN_DIAGNOSTIC_COMPLETE")
