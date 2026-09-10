options(stringsAsFactors=FALSE,width=150);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-"/root/projects/UC_Treatment_Recovery";setwd(root)
out<-"runs/001C2_refinement"
o<-readRDS("runs/001C_programs/program_scores.rds");p<-o$pairs;D<-o$delta;BL<-o$baseline
focus<-c("MATURE_COLONOCYTE6_EXPLORATORY","UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy","UC_IA13_JAKSTAT_2024_gene_mean_proxy","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
res<-list()
for(s in focus){
 dd<-p;dd$delta<-D[s,];dd$baseline<-BL[s,];dd$inflam<-D["HALLMARK_INFLAMMATORY_RESPONSE",]
 for(spec in c("common_baseline","common_baseline_inflammation")){
  f<-lm(if(spec=="common_baseline")delta~arm+response+baseline else delta~arm+response+baseline+inflam,data=dd)
  X<-model.matrix(f);e<-residuals(f);h<-hatvalues(f);bread<-solve(crossprod(X))
  # HC3 accounts for leverage and unequal residual variances; t critical value uses residual df and remains approximate.
  meat<-crossprod(X*as.numeric(e/(1-h)))
  V<-bread%*%meat%*%bread
  j<-match("responseYes",colnames(X));se<-sqrt(V[j,j]);est<-unname(coef(f)[j]);df<-df.residual(f)
  Xs<-cbind(Intercept=1,scale(X[,-1,drop=FALSE]))
  residual_sd_ratio<-with(dd,max(tapply(e,response,sd))/min(tapply(e,response,sd)))
  res[[length(res)+1]]<-data.table(program=s,model=spec,estimate=est,HC3_se=se,HC3_lower=est-qt(.975,df)*se,HC3_upper=est+qt(.975,df)*se,HC3_p=2*pt(-abs(est/se),df),df=df,scaled_design_condition_number=kappa(Xs,exact=TRUE),max_leverage=max(h),residual_response_sd_ratio=residual_sd_ratio)
 }
}
rr<-rbindlist(res);rr[,HC3_FDR:=p.adjust(HC3_p,"BH"),by=model]
fwrite(rr,file.path(out,"HC3_common_model_sensitivity.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(out,"HC3_sessionInfo.txt"))
cat("HC3_COMPLETE\n")
