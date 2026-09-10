options(stringsAsFactors=FALSE,width=170);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(limma);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/002A_GSE23597";dir.create(out,showWarnings=FALSE,recursive=TRUE)
lock<-read_json("signature_lock.json",simplifyVector=TRUE)
stopifnot(lock$version=="1.1",lock$status=="FROZEN")
fail<-function(reason){write_json(list(status="COMPLETE",primary_evaluable=FALSE,reason=reason,cohort="GSE23597",alternative_primary_allowed=FALSE),file.path(out,"validation_summary.json"),auto_unbox=TRUE,pretty=TRUE);stop(reason)}
con<-gzfile("inputs/validation/GSE23597_series_matrix.txt.gz","rt")
repeat{line<-readLines(con,n=1,warn=FALSE);if(!length(line))stop("No matrix start");if(line=="!series_matrix_table_begin")break}
x<-as.matrix(read.delim(con,check.names=FALSE,quote='"',comment.char="!",row.names=1));close(con);storage.mode(x)<-"double"
if(any(!is.finite(x)))fail("Nonfinite input matrix; source clarification needed")
if(anyDuplicated(rownames(x))||anyDuplicated(colnames(x)))fail("Duplicate source matrix row/sample IDs")
q<-quantile(x,c(0,.01,.25,.5,.75,.99,1))
floored<-0L;transform<-"retained deposited log scale"
if(q[6]>100&&q[1]>=0){floored<-sum(x<1);x<-log2(pmax(x,1));transform<-"GCOS globally scaled intensity -> log2(pmax(x,1)); no second normalization"} else if(!(q[6]<=40&&q[1]>=-20))fail("Scale outside predeclared technical rules")
registry<-fread("runs/001A_metadata/sample_registry_v1.tsv",na.strings="")
m<-registry[series=="GSE23597"];m<-m[match(colnames(x),gsm)]
stopifnot(identical(m$gsm,colnames(x)))
pairfile<-"provenance/reviews/002A_metadata_prevalidation_pairs_independent_v1.tsv"
p<-as.data.frame(fread(pairfile));stopifnot(nrow(p)==32,!anyDuplicated(p$subject),!"P13"%in%p$subject,all(p$baseline_gsm%in%colnames(x)),all(p$week8_gsm%in%colnames(x)))
p$response<-factor(p$response,levels=c("No","Yes"));p$treatment_dose<-factor(p$treatment_dose,levels=c("Placebo","IFX_5mgkg","IFX_10mgkg"))
map<-fread("runs/001D_annotation/GPL570_canonical_probe_mapping.tsv",na.strings="")
map<-map[match(rownames(x),probe)]
keep<-!is.na(map$retained)&map$retained==1&!is.na(map$canonical_symbol)
g<-avereps(x[keep,,drop=FALSE],ID=map$canonical_symbol[keep]);g<-g[order(rownames(g)),,drop=FALSE]
cand<-lock$candidate$id;infl<-"HALLMARK_INFLAMMATORY_RESPONSE";common<-lock$inflammation_covariate$genes
cgenes<-lock$candidate$genes
mem<-fread("runs/001D_annotation/program_membership_canonical.tsv")
refs<-lock$secondary_reference_programs
stopifnot(setequal(mem[program==cand]$gene,cgenes),length(common)==194L,!anyDuplicated(common))
selected<-c(cand,infl,refs)
cov<-list();sets<-list()
for(s in selected){
 S<-unique(mem[program==s]$gene);hits<-intersect(S,rownames(g));sets[[s]]<-hits
 cov[[length(cov)+1]]<-data.table(program=s,source_n=length(S),measured_n=length(hits),coverage=length(hits)/length(S),missing=paste(setdiff(S,hits),collapse=";"))
}
fwrite(rbindlist(cov),file.path(out,"actual_program_coverage.tsv"),sep="\t")
if(!all(cgenes%in%rownames(g)))fail("Candidate6 incomplete")
if(!all(common%in%rownames(g)))fail("Frozen common194 inflammation panel incomplete")
if(mean(mem[program==infl]$gene%in%rownames(g))<.8)fail("Original inflammation200 coverage below80%")
sets[[infl]]<-common
fwrite(data.table(gene=common,measured=common%in%rownames(g),weight=1),file.path(out,"actual_scored_inflammation194.tsv"),sep="\t")
S<-t(vapply(sets,function(z)colMeans(g[z,,drop=FALSE]),numeric(ncol(g))))
colnames(S)<-colnames(g)
BL<-S[,p$baseline_gsm];PO<-S[,p$week8_gsm];D<-PO-BL;colnames(BL)<-colnames(PO)<-colnames(D)<-p$subject
saveRDS(list(probe_expression_log2=x,gene_expression=g,metadata=m,pairs=p,sample_scores=S,baseline=BL,post=PO,delta=D,sets=sets),file.path(out,"validation_scores.rds"))
dd<-p;dd$delta_candidate<-D[cand,];dd$baseline_candidate<-BL[cand,];dd$delta_inflammatory_response<-D[infl,]
fwrite(dd,file.path(out,"primary_patient_data.tsv"),sep="\t")
fullform<-as.formula(lock$primary_validation$model)
mf<-model.frame(fullform,dd,na.action=na.pass)
if(any(!complete.cases(mf)))fail("Essential fixed model field missing")
fit<-lm(fullform,dd)
if(fit$rank!=length(coef(fit))||!"responseYes"%in%names(coef(fit))||any(!is.finite(coef(fit))))fail("Primary design rank deficient or response not estimable")
rowfit<-function(fit,s,model){
 a<-summary(fit)$coef["responseYes",];ci<-confint(fit,"responseYes");X<-model.matrix(fit);Z<-cbind(1,scale(X[,-1,drop=FALSE]));vv<-diag(solve(cor(X[,-1,drop=FALSE])))
 data.table(program=s,model=model,estimate=unname(a[1]),se=unname(a[2]),lower=ci[1],upper=ci[2],p=unname(a[4]),n=nobs(fit),df=fit$df.residual,rank=fit$rank,raw_condition=kappa(X,exact=TRUE),standardized_condition=kappa(Z,exact=TRUE),max_VIF=max(vv),max_leverage=max(hatvalues(fit)))
}
hc3<-function(fit,s,model){
 X<-model.matrix(fit);h<-hatvalues(fit);e<-resid(fit);bread<-solve(crossprod(X));v<-bread%*%crossprod(X,X*as.numeric((e/(1-h))^2))%*%bread;be<-coef(fit)["responseYes"];se<-sqrt(v["responseYes","responseYes"]);df<-df.residual(fit)
 data.table(program=s,model=model,estimate=unname(be),se=unname(se),lower=unname(be-qt(.975,df)*se),upper=unname(be+qt(.975,df)*se),p=unname(2*pt(-abs(be/se),df)),n=nobs(fit),df=df)
}
primary<-rowfit(fit,cand,"primary_baseline_inflammation")
fwrite(primary,file.path(out,"primary_test.tsv"),sep="\t")
models<-list(primary_baseline_inflammation=fit,mandatory_no_baseline=lm(as.formula(lock$primary_validation$mandatory_parallel_model$model),dd),without_inflammation=lm(delta_candidate~response+treatment_dose+baseline_candidate,dd))
sub<-droplevels(dd[dd$treatment_dose!="Placebo",]);ifit<-lm(fullform,sub)
if(ifit$rank==length(coef(ifit)))models$IFX_only<-ifit
effects<-rbindlist(lapply(names(models),function(z)rowfit(models[[z]],cand,z)))
rob<-rbindlist(lapply(names(models),function(z)hc3(models[[z]],cand,z)))
fwrite(effects,file.path(out,"candidate_models.tsv"),sep="\t")
fwrite(rob,file.path(out,"candidate_models_HC3.tsv"),sep="\t")
diag<-data.table(subject=dd$subject,response=dd$response,treatment_dose=dd$treatment_dose,leverage=hatvalues(fit),residual=resid(fit),studentized_residual=rstudent(fit),cooks_distance=cooks.distance(fit))
fwrite(diag,file.path(out,"primary_patient_diagnostics.tsv"),sep="\t")
X<-model.matrix(fit)[,-1,drop=FALSE];fwrite(data.table(term=colnames(X),VIF=diag(solve(cor(X)))),file.path(out,"primary_VIF.tsv"),sep="\t")
lopo<-rbindlist(lapply(seq_len(nrow(dd)),function(i){f<-lm(fullform,dd[-i,]);data.table(omitted_subject=dd$subject[i],estimate=unname(coef(f)["responseYes"]),p=summary(f)$coef["responseYes",4])}))
fwrite(lopo,file.path(out,"leave_one_patient_out.tsv"),sep="\t")
secondary<-list();sec_no_baseline<-list()
for(s in refs){
 cc<-rbindlist(cov)[program==s]
 if(cc$coverage<.8){secondary[[length(secondary)+1]]<-data.table(program=s,model="secondary_baseline_inflammation",evaluable=FALSE,p=NA_real_);next}
 rr<-dd;rr$delta_candidate<-D[s,];rr$baseline_candidate<-BL[s,]
 f<-lm(fullform,rr)
 if(f$rank!=length(coef(f))||any(!is.finite(coef(f)))){secondary[[length(secondary)+1]]<-data.table(program=s,model="secondary_baseline_inflammation",evaluable=FALSE,p=NA_real_,reason="secondary model rank deficient");next}
 a<-rowfit(f,s,"secondary_baseline_inflammation");a$evaluable<-TRUE;secondary[[length(secondary)+1]]<-a
 sec_no_baseline[[length(sec_no_baseline)+1]]<-rowfit(lm(as.formula(lock$primary_validation$mandatory_parallel_model$model),rr),s,"reference_no_baseline_descriptive")
}
ss<-rbindlist(secondary,fill=TRUE);ss$FDR_across_three<-p.adjust(ss$p,method="BH",n=3)
fwrite(ss,file.path(out,"secondary_reference_tests.tsv"),sep="\t")
fwrite(rbindlist(sec_no_baseline),file.path(out,"secondary_no_baseline_descriptive.tsv"),sep="\t")
saveRDS(models,file.path(out,"candidate_model_fits.rds"))
qc<-data.table(gsm=colnames(x),median_log2=apply(x,2,median),min_log2=apply(x,2,min),max_log2=apply(x,2,max),sd_log2=apply(x,2,sd))
fwrite(qc,file.path(out,"sample_qc.tsv"),sep="\t")
dups<-colnames(x)[duplicated(as.data.frame(t(x)))]
qcsummary<-list(probes=nrow(x),arrays=ncol(x),canonical_genes=nrow(g),source_quantiles=as.list(q),log_transform=transform,floored_values=floored,exact_duplicate_expression_arrays=dups,source_metadata_matched=TRUE,pairs=nrow(p),patient_exclusion_by_expression=FALSE,all_common194_present=TRUE)
write_json(qcsummary,file.path(out,"qc_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=15)
info<-list(status="COMPLETE",cohort="GSE23597",primary_evaluable=TRUE,primary=as.list(primary),primary_support_rule_met=primary$estimate>0&&primary$p<.05,primary_tests=1,secondary_tests=3,discovery_genes_unchanged=TRUE,validation_expression_accessed=TRUE,alternative_primary_allowed=FALSE,claim="Concurrent conditional association only; no baseline prediction, causal repair, novel cell state or comparator superiority",lock_version=lock$version,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
write_json(info,file.path(out,"validation_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=15)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(effects);print(rob);print(ss);cat("002A_VALIDATION_COMPLETE\n")
