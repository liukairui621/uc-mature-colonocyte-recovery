options(stringsAsFactors=FALSE,width=160);set.seed(20260911)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery")
out<-"runs/002D_VDZ_support";dir.create(out,showWarnings=FALSE,recursive=TRUE)
if(file.exists(file.path(out,"candidate_models.tsv"))||file.exists(file.path(out,"analysis_summary.json")))stop("Existing002Dresults: refuse overwrite")
plan<-read_json("planning/support_plan_002D_VDZ.frozen.json",simplifyVector=TRUE)
stopifnot(plan$status=="FROZEN",!plan$can_replace_GSE23597)
for(i in seq_len(nrow(plan$inputs))){
 f<-plan$inputs$path[i];h<-strsplit(system2("sha256sum",f,stdout=TRUE)," ")[[1]][1]
 stopifnot(identical(h,plan$inputs$sha256[i]))
}
z<-readRDS("runs/002B_GSE73661_support/support_scores.rds")
lock<-read_json("signature_lock.json",simplifyVector=TRUE)
cand<-lock$candidate$id;infl<-"HALLMARK_INFLAMMATORY_RESPONSE";refs<-lock$secondary_reference_programs
stopifnot(setequal(z$sets[[cand]],lock$candidate$genes),setequal(z$sets[[infl]],lock$inflammation_covariate$genes),length(z$sets[[infl]])==194)
p<-fread("runs/002D_VDZ_metadata/VDZ_pairs.tsv")
stopifnot(nrow(p)==40,!anyDuplicated(p$subject),all(c(p$baseline_gsm,p$followup_gsm)%in%colnames(z$sample_scores)))
p$endoscopic_healing<-factor(p$endoscopic_healing,levels=c("No","Yes"))
allrows<-list();diags<-list();loo<-list();patients<-list();refrows<-list();fitstore<-list()
fitrow<-function(f,visit,program,model,method="OLS"){
 X<-model.matrix(f);term<-"endoscopic_healingYes";N<-nobs(f);df<-df.residual(f)
 ok<-f$rank==ncol(X)&&df>=5&&term%in%rownames(summary(f)$coefficients)
 if(!ok)return(data.table(visit=visit,program=program,model=model,method=method,evaluable=FALSE,n=N,df=df,reason="Rank deficiency, term missing, or residual df below5",p=NA_real_))
 b<-unname(coef(f)[term]);se<-summary(f)$coef[term,2]
 if(method=="HC3"){
  h<-hatvalues(f);if(any(h>=1-1e-10))return(data.table(visit=visit,program=program,model=model,method=method,evaluable=FALSE,n=N,df=df,reason="HC3 leverage one",p=NA_real_))
  bread<-solve(crossprod(X));vv<-bread%*%crossprod(X,X*as.numeric((resid(f)/(1-h))^2))%*%bread
  se<-sqrt(vv[term,term])
 }
 data.table(visit=visit,program=program,model=model,method=method,evaluable=TRUE,n=N,df=df,estimate=b,se=se,lower=b-qt(.975,df)*se,upper=b+qt(.975,df)*se,p=2*pt(-abs(b/se),df),max_leverage=max(hatvalues(f)),residual_sigma=summary(f)$sigma,reason="")
}
for(visit_i in c("W6","W12")){
 d<-as.data.frame(p[visit==visit_i])
 stopifnot(min(table(d$endoscopic_healing))>=2)
 BL<-z$sample_scores[,d$baseline_gsm,drop=FALSE];PO<-z$sample_scores[,d$followup_gsm,drop=FALSE];D<-PO-BL
 d$baseline_candidate<-as.numeric(BL[cand,]);d$post_candidate<-as.numeric(PO[cand,]);d$delta_candidate<-as.numeric(D[cand,]);d$delta_common194_inflammation<-as.numeric(D[infl,])
 patients[[visit_i]]<-as.data.table(d)
 fits<-list()
 for(label in names(plan$models)){
  f<-lm(as.formula(plan$models[[label]]),d,na.action=na.fail);fits[[label]]<-f
  for(meth in c("OLS","HC3"))allrows[[length(allrows)+1]]<-fitrow(f,visit_i,cand,label,meth)
 }
 fitstore[[visit_i]]<-fits
 f<-fits$adjusted_support
 diags[[visit_i]]<-data.table(visit=visit_i,subject=d$subject,endoscopic_healing=d$endoscopic_healing,leverage=hatvalues(f),cooks_distance=cooks.distance(f),residual=resid(f))
 loo[[visit_i]]<-rbindlist(lapply(seq_len(nrow(d)),function(i){
  f1<-lm(formula(f),d[-i,],na.action=na.fail);rr<-fitrow(f1,visit_i,cand,"leave_one_patient_out");rr$omitted_subject<-d$subject[i];rr
 }),fill=TRUE)
 for(s in refs){
  dr<-d;dr$baseline_candidate<-as.numeric(BL[s,]);dr$delta_candidate<-as.numeric(D[s,])
  fr<-lm(as.formula(plan$models$adjusted_support),dr,na.action=na.fail)
  refrows[[length(refrows)+1]]<-fitrow(fr,visit_i,s,"reference_descriptive")
 }
}
effects<-rbindlist(allrows,fill=TRUE);refeffects<-rbindlist(refrows,fill=TRUE)
refeffects$BH_across_six_descriptive<-p.adjust(refeffects$p,"BH",n=6)
fwrite(effects,file.path(out,"candidate_models.tsv"),sep="\t")
fwrite(refeffects,file.path(out,"reference_models.tsv"),sep="\t")
fwrite(rbindlist(patients),file.path(out,"patient_level_scores.tsv"),sep="\t")
fwrite(rbindlist(diags),file.path(out,"patient_diagnostics.tsv"),sep="\t")
fwrite(rbindlist(loo,fill=TRUE),file.path(out,"leave_one_patient_out.tsv"),sep="\t")
saveRDS(fitstore,file.path(out,"model_fits.rds"))
main<-effects[visit=="W6"&model=="adjusted_support"&method=="OLS"]
sec<-effects[visit=="W12"&model=="adjusted_support"&method=="OLS"]
classification<-if(!main$evaluable)"not_evaluable" else if(main$estimate>0&&main$p<.05)"positive_supportive_association" else "cross_drug_support_not_established"
su<-list(status="COMPLETE",role=plan$role,main_W6=as.list(main),secondary_W12=as.list(sec),interpretation=classification,previous_VDZ_sample_scores_computed=TRUE,original_GSE23597_primary_support=FALSE,pooled_P_computed=FALSE,drug_comparison_computed=FALSE,negative_result_does_not_prove_anti_TNF_specificity=TRUE,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
write_json(su,file.path(out,"analysis_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(effects);print(refeffects);cat("002D_VDZ_COMPLETE\n")
