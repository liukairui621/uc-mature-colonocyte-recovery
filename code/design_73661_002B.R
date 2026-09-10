options(stringsAsFactors=FALSE);suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/002B_precision"
p<-fread(file.path(out,"GSE73661_IFX_pairs_metadata_independent_v1.tsv"))
stopifnot(nrow(p)==23,all(p$treatment_path=="IFX"),all(p$followup_visit=="W4_W6"),sum(p$endoscopic_healing=="Yes")==8)
N<-nrow(p);prop<-mean(p$endoscopic_healing=="Yes");df<-N-4
fitted<-fread(file.path(out,"fitted_precision_components.tsv"))
power_t<-function(effect,se,df){z<-qt(.975,df);pt(-z,df,ncp=effect/se)+pt(z,df,ncp=effect/se,lower.tail=FALSE)}
mde<-function(se,df)uniroot(function(e)power_t(e,se,df)-.8,c(0,20*se),tol=1e-11)$root
grid<-CJ(sigma_source=c("GSE92415","GSE23597"),assumed_response_VIF=c(1,1.25,1.5,2))
grid$residual_sigma<-fitted$residual_sigma[match(grid$sigma_source,fitted$cohort)]
grid$SE<-grid$residual_sigma*sqrt(grid$assumed_response_VIF/(N*prop*(1-prop)))
grid$MDE80<-vapply(grid$SE,mde,numeric(1),df=df)
grid$power_at_discovery_0_599<-vapply(grid$SE,power_t,numeric(1),effect=fitted[cohort=="GSE92415"]$estimate,df=df)
grid$n<-N;grid$n_healing<-8;grid$n_nonhealing<-15;grid$expected_df<-df
fwrite(grid,file.path(out,"GSE73661_preexpression_precision_scenarios.tsv"),sep="\t")
# Replicate current response proportion and VIF only to illustrate planning; do not claim true future precision.
ngrid<-CJ(n=c(23,40,60,80,100,150,200),sigma_source=c("GSE92415","GSE23597"),assumed_response_VIF=c(1.25,1.5))
ngrid$residual_sigma<-fitted$residual_sigma[match(ngrid$sigma_source,fitted$cohort)]
ngrid$assumed_response_fraction<-prop
ngrid$SE<-ngrid$residual_sigma*sqrt(ngrid$assumed_response_VIF/(ngrid$n*prop*(1-prop)))
ngrid$MDE80<-mapply(mde,ngrid$SE,ngrid$n-4)
ngrid$power_at_0_599<-mapply(power_t,effect=fitted[cohort=="GSE92415"]$estimate,se=ngrid$SE,df=ngrid$n-4)
fwrite(ngrid,file.path(out,"future_sample_size_scenarios.tsv"),sep="\t")
write_json(list(status="COMPLETE",GSE73661_expression_accessed=FALSE,phase="After GSE23597 results but before GSE73661 expression",cohort="GSE73661 IFX only",n=23,endoscopic_healing=8,nonhealing=15,model="delta_candidate ~ endoscopic_healing + baseline_candidate + delta_common194_inflammation",expected_parameters=4,expected_residual_df=19,actual_full_model_covariance_unknown=TRUE,role="Descriptive/supportive-only regardless of these scenarios or later results; never replaces GSE23597",assumptions=c("Transport residual SD values as illustrative scenarios only","Use VIF grid1,1.25,1.5,2 because actual expression covariates unseen","No claim GPL6244 residual noise equals RMA discovery or that array type identifies variance","Clinical response and endoscopic healing are distinct endpoints; effect transport0.599 is an illustrative assumption","Visit and numeric dose unavailable; do not invent them"),completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE)),file.path(out,"GSE73661_preexpression_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=15)
print(grid);cat("GSE73661_PRECISION_COMPLETE\n")
