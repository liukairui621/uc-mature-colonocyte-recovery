options(stringsAsFactors=FALSE,width=160);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-"/root/projects/UC_Treatment_Recovery";setwd(root)
out<-"runs/001C2_refinement";dir.create(out,showWarnings=FALSE,recursive=TRUE)
o<-readRDS("runs/001C_programs/program_scores.rds");p<-o$pairs;D<-o$delta;BL<-o$baseline
candidate<-"MATURE_COLONOCYTE6_EXPLORATORY"
CT<-"UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy"
infl<-"HALLMARK_INFLAMMATORY_RESPONSE";epi<-"EPITHELIAL_ABUNDANCE_PROXY4"
focus<-c(candidate,CT,"UC_IA13_JAKSTAT_2024_gene_mean_proxy","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
p$arm<-factor(p$arm,levels=c("Placebo","golimumab"));p$response<-factor(p$response,levels=c("No","Yes"))
batch<-fread("provenance/reviews/metadata_independent_v1_gsm_barcode.tsv",sep="\t")
bb<-batch$inferred_cel_barcode[match(p$baseline_gsm,batch$gsm)];bp<-batch$inferred_cel_barcode[match(p$post_gsm,batch$gsm)]
lv<-sort(unique(c(bb,bp)));B<-sapply(lv[-1],function(z)as.integer(bp==z)-as.integer(bb==z));colnames(B)<-paste0("barcode",seq_len(ncol(B)))
res<-list();lopo<-list();logistic<-list()
cirow<-function(fit,program,model){
 a<-summary(fit)$coef["responseYes",];ci<-confint(fit,"responseYes")
 data.table(program=program,model=model,estimate=unname(a[1]),se=unname(a[2]),lower=ci[1],upper=ci[2],p=unname(a[4]),n=nobs(fit),df=fit$df.residual,rank=fit$rank,design_condition_number=kappa(model.matrix(fit),exact=TRUE))
}
for(s in focus){
 dd<-p;dd$delta<-D[s,];dd$baseline<-BL[s,];dd$inflam<-D[infl,];dd$epi<-D[epi,];dd$ct_delta<-D[CT,];dd$ct_baseline<-BL[CT,]
 dd<-cbind(dd,B)
 common<-c("arm","response","baseline")
 models<-list(common_baseline=common,common_baseline_inflammation=c(common,"inflam"),common_baseline_inflammation_epithelial=c(common,"inflam","epi"),common_baseline_inflammation_barcode=c(common,"inflam",colnames(B)))
 if(s%in%c(candidate,CT)){
  for(src in c("MARS_TableS1_column_","MARS_MSigDB74_")){
   sn<-grep(paste0("^",src),rownames(D),value=TRUE)
   mn<-paste0(if(src=="MARS_TableS1_column_")"marsT" else "marsV",seq_along(sn))
   for(j in seq_along(mn))dd[[mn[j]]]<-D[sn[j],]
   models[[paste0("common_baseline_inflammation_",if(src=="MARS_TableS1_column_")"MARS_TableS1_11" else "MARS_v74_11")]]<-c(common,"inflam",mn)
  }
 }
 if(s==candidate)models$candidate_conditioned_on_CT_program<-c(common,"inflam","ct_delta","ct_baseline")
 for(model in names(models)){
  f<-lm(reformulate(models[[model]],"delta"),data=dd)
  stopifnot(f$rank==length(coef(f)),all(is.finite(coef(f))))
  res[[length(res)+1]]<-cirow(f,s,model)
  if(model=="common_baseline_inflammation")for(i in seq_len(nrow(dd))){
   ff<-lm(reformulate(models[[model]],"delta"),data=dd[-i,])
   lopo[[length(lopo)+1]]<-data.table(program=s,omitted_subject=dd$subject[i],estimate=unname(coef(ff)["responseYes"]))
  }
 }
 # Adjusted contemporaneous association, no prediction or model selection claim.
 dd$resp<-as.integer(dd$response=="Yes")
 f0<-glm(resp~arm+baseline+inflam,family=binomial(),data=dd)
 f1<-glm(resp~arm+baseline+inflam+delta,family=binomial(),data=dd)
 a<-summary(f1)$coef["delta",];lrt<-anova(f0,f1,test="LRT");sdx<-sd(dd$delta);be<-unname(a[1]);se<-unname(a[2])
 logistic[[length(logistic)+1]]<-data.table(program=s,beta=be,se=se,OR_per_SD=exp(be*sdx),lower_OR=exp((be-qnorm(.975)*se)*sdx),upper_OR=exp((be+qnorm(.975)*se)*sdx),LRT_p=lrt$"Pr(>Chi)"[2],converged=f1$converged)
}
ee<-rbindlist(res);ee[,FDR_across_declared_programs:=p.adjust(p,"BH"),by=model]
fwrite(ee,file.path(out,"refinement_response_effects_CI.tsv"),sep="\t")
fwrite(rbindlist(lopo),file.path(out,"refinement_leave_one_patient_out.tsv"),sep="\t")
ll<-rbindlist(logistic);ll[,FDR_across_four_programs:=p.adjust(LRT_p,"BH")]
fwrite(ll,file.path(out,"common_covariate_incremental_association.tsv"),sep="\t")
# Quantify score-level technical proxy sensitivity without interpreting barcode as verified batch.
armrows<-list()
for(s in c(candidate,"HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_OXIDATIVE_PHOSPHORYLATION")){
 dd<-cbind(p,B);dd$delta<-D[s,]
 for(branch in c("unadjusted","barcode_proxy")){
  fm<-if(branch=="unadjusted")lm(delta~0+arm,data=dd) else lm(reformulate(c("arm",colnames(B)),"delta",intercept=FALSE),data=dd)
  co<-summary(fm)$coef;ci<-confint(fm)
  for(term in paste0("arm",levels(p$arm)))armrows[[length(armrows)+1]]<-data.table(program=s,branch=branch,arm=sub("^arm","",term),estimate=co[term,1],se=co[term,2],lower=ci[term,1],upper=ci[term,2],p=co[term,4])
 }
}
fwrite(rbindlist(armrows),file.path(out,"barcode_program_change_sensitivity.tsv"),sep="\t")
write_json(list(status="COMPLETE",scope="Post-observation refinement; compare existing mature-cell program and review-proposed six-gene candidate; no new-state claim",hypothesis_timing="Added after observing candidate correlation with CT program in 001C; all results exploratory",common_covariates="Treatment arm and own baseline score; inflammation-score change in conditional association; baseline total Mayo omitted here because unavailable in GSE23597, retained in 001C sensitivity",MARS_limitation="11 gene-mean component proxies; not original GSVA/clusters",IBrD_excluded_from_formal_adjustment="Only15/81genes; inadequate coverage and pretrained model unavailable",candidate_frozen=FALSE,validation_expression_accessed=FALSE),file.path(out,"refinement_summary.json"),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat("REFINEMENT_COMPLETE\n")
