options(stringsAsFactors=FALSE,width=160);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/002B_precision";dir.create(out,showWarnings=FALSE,recursive=TRUE)
cand<-"MATURE_COLONOCYTE6_EXPLORATORY"
d<-readRDS("runs/001D_annotation/discovery_model_fits.rds")[[paste(cand,"primary_baseline_inflammation",sep="::")]]
v<-readRDS("runs/002A_GSE23597/candidate_model_fits.rds")$primary_baseline_inflammation
power_t<-function(effect,se,df,alpha=.05){crit<-qt(1-alpha/2,df);pt(-crit,df,ncp=effect/se)+pt(crit,df,ncp=effect/se,lower.tail=FALSE)}
mde<-function(se,df,target=.8)uniroot(function(delta)power_t(delta,se,df)-target,c(0,se*20),tol=1e-11)$root
model_info<-function(f,cohort){
 X<-model.matrix(f);j<-match("responseYes",colnames(X));r<-X[,j];Z<-X[,-j,drop=FALSE]
 residual_response<-lm.fit(Z,r)$residuals;SST<-sum((r-mean(r))^2);RSSr<-sum(residual_response^2)
 sigma<-summary(f)$sigma;VIF<-SST/RSSr;se<-summary(f)$coef[j,2];df<-df.residual(f)
 stopifnot(abs(se-sigma/sqrt(RSSr))<1e-10)
 data.table(cohort=cohort,n=nobs(f),p_response=mean(r),df=df,delta_SD=sd(model.response(model.frame(f))),residual_sigma=sigma,residual_sigma_lower=sqrt(df*sigma^2/qchisq(.975,df)),residual_sigma_upper=sqrt(df*sigma^2/qchisq(.025,df)),response_SST=SST,residualized_response_SS=RSSr,response_VIF=VIF,SE=se,estimate=unname(coef(f)[j]),MDE80_plugin=mde(se,df),MDE80_lower_sigma=mde(se*sqrt(df/qchisq(.975,df)),df),MDE80_upper_sigma=mde(se*sqrt(df/qchisq(.025,df)),df))
}
info<-rbindlist(list(model_info(d,"GSE92415"),model_info(v,"GSE23597")))
fwrite(info,file.path(out,"fitted_precision_components.tsv"),sep="\t")
di<-info[cohort=="GSE92415"];vi<-info[cohort=="GSE23597"]
factors<-data.table(component=c("sample_size","response_balance","response_collinearity","residual_sigma"),factor=c(sqrt(di$n/vi$n),sqrt(di$p_response*(1-di$p_response)/(vi$p_response*(1-vi$p_response))),sqrt(vi$response_VIF/di$response_VIF),vi$residual_sigma/di$residual_sigma))
prod_factors<-prod(factors$factor);actual<-vi$SE/di$SE
stopifnot(abs(prod_factors-actual)<1e-10)
fwrite(factors,file.path(out,"exact_SE_ratio_decomposition.tsv"),sep="\t")
curve<-rbindlist(lapply(seq_len(nrow(info)),function(i){r<-info[i];data.table(cohort=r$cohort,assumed_true_effect=seq(0,1.8,by=.01),power=vapply(seq(0,1.8,by=.01),power_t,numeric(1),se=r$SE,df=r$df),SE_plugin=r$SE,df=r$df)}))
fwrite(curve,file.path(out,"conditional_power_scenarios.tsv"),sep="\t")
effect_scenarios<-c(.25,.4,di$estimate,.8,1)
powerrows<-rbindlist(lapply(effect_scenarios,function(e)data.table(assumed_effect=e,validation_power=power_t(e,vi$SE,vi$df),label=if(e==di$estimate)"Discovery-selected point estimate; potentially optimistic" else "Illustrative fixed effect assumption")))
fwrite(powerrows,file.path(out,"selected_power_scenarios.tsv"),sep="\t")
observed<-data.table(effect=vi$estimate,observed_power=power_t(vi$estimate,vi$SE,vi$df),role="Audit of reviewer arithmetic only; not independent evidence of why P was nonsignificant")
fwrite(observed,file.path(out,"observed_power_arithmetic_only.tsv"),sep="\t")
# These are biological-plus-technical total dispersions, not isolated measurement error.
od<-readRDS("runs/001D_annotation/discovery_preprocessed_canonical.rds")
ov<-readRDS("runs/002A_GSE23597/validation_scores.rds")
genes<-c("AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12")
gd<-od$delta[genes,,drop=FALSE]
gv<-ov$gene_expression[genes,ov$pairs$week8_gsm,drop=FALSE]-ov$gene_expression[genes,ov$pairs$baseline_gsm,drop=FALSE]
sdtab<-data.table(feature=c(genes,"Six_gene_mean"),discovery_delta_SD=c(apply(gd,1,sd),sd(colMeans(gd))),validation_delta_SD=c(apply(gv,1,sd),sd(colMeans(gv))))
sdtab$ratio<-sdtab$validation_delta_SD/sdtab$discovery_delta_SD
fwrite(sdtab,file.path(out,"observed_delta_dispersion.tsv"),sep="\t")
# Variance of a mean contains covariances; average scores do not prove amplified platform noise.
covrows<-list()
for(co in c("discovery","validation")){
 G<-if(co=="discovery")gd else gv;C<-cov(t(G));k<-nrow(C)
 covrows[[co]]<-data.table(cohort=co,variance_of_score_mean=var(colMeans(G)),sum_individual_variances_over_k2=sum(diag(C))/k^2,sum_cross_covariances_over_k2=(sum(C)-sum(diag(C)))/k^2,mean_gene_pair_correlation=mean(cor(t(G))[upper.tri(C)]))
}
fwrite(rbindlist(covrows),file.path(out,"score_variance_identity.tsv"),sep="\t")
# Descriptive coefficient attenuation in distinct observational cohorts; no mediated proportion.
dr<-readRDS("runs/001D_annotation/discovery_model_fits.rds")
vr<-readRDS("runs/002A_GSE23597/candidate_model_fits.rds")
at<-data.table(cohort=c("GSE92415","GSE23597"),before=c(unname(coef(dr[[paste(cand,"baseline_only",sep="::")]])["responseYes"]),unname(coef(vr$without_inflammation)["responseYes"])),after=c(di$estimate,vi$estimate))
at$attenuation_percent<-100*(1-at$after/at$before)
fwrite(at,file.path(out,"inflammation_conditioning_attenuation.tsv"),sep="\t")
loo<-fread("runs/002A_GSE23597/leave_one_patient_out.tsv")
patient<-fread("runs/002A_GSE23597/primary_patient_diagnostics.tsv");imax<-which.max(patient$cooks_distance)
summary<-list(status="COMPLETE",statistical_role="Post-validation precision audit and conditional design sensitivity, not post-hoc proof of why the result was negative",validation_primary_unchanged=TRUE,validation_MDE80_plugin=vi$MDE80_plugin,power_at_discovery_selected_estimate=power_t(di$estimate,vi$SE,vi$df),observed_power_not_for_inference=observed$observed_power,actual_SE_ratio=actual,exact_factor_product=prod_factors,LOPO=list(min_estimate=min(loo$estimate),max_estimate=max(loo$estimate),min_p=min(loo$p),max_p=max(loo$p)),max_Cooks_D=list(subject=patient$subject[imax],value=patient$cooks_distance[imax]),assumptions=c("Fixed observed covariate design and plug-in residual SD under homoscedastic normal linear model","Power at discovery-selected estimate may be optimistic because of selection and cross-cohort effect differences","Residual-sigma intervals condition on model assumptions; not cross-platform causal evidence","Raw delta SD combines biology,treatment,time,patient mix,technical processing and random variation","No paired same biospecimens processed both ways; cannot identify RMA->GCOS noise multiplier","No new pooled P value or rescue primary validation"),completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
write_json(summary,file.path(out,"precision_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=15)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(info);print(factors);print(powerrows);print(sdtab);print(at);cat("PRECISION_COMPLETE\n")
