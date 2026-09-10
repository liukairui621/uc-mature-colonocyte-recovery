options(stringsAsFactors=FALSE,width=160);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(limma);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/002B_GSE73661_support";dir.create(out,recursive=TRUE,showWarnings=FALSE)
plan<-read_json("planning/support_plan_002B.frozen.json",simplifyVector=TRUE)
stopifnot(plan$status=="FROZEN",plan$role=="descriptive_supportive_only",!plan$can_replace_GSE23597)
lock<-read_json("signature_lock.json",simplifyVector=TRUE)
fail<-function(reason){write_json(list(status="COMPLETE",evaluable=FALSE,reason=reason,role=plan$role,can_replace_GSE23597=FALSE),file.path(out,"support_summary.json"),pretty=TRUE,auto_unbox=TRUE);stop(reason)}
con<-gzfile("inputs/support/GSE73661_series_matrix.txt.gz","rt");headers<-character()
repeat{line<-readLines(con,n=1,warn=FALSE);if(!length(line))fail("Matrix start absent");if(line=="!series_matrix_table_begin")break;headers<-c(headers,line)}
x<-as.matrix(read.delim(con,check.names=FALSE,quote='"',comment.char="!",row.names=1));close(con);storage.mode(x)<-"double"
writeLines(headers,file.path(out,"downloaded_matrix_header.txt"))
if(any(!is.finite(x))||anyDuplicated(colnames(x))||anyDuplicated(rownames(x)))fail("Nonfinite or duplicate matrix IDs")
q<-quantile(x,c(0,.01,.25,.5,.75,.99,1))
if(q[6]>40||q[1]< -20)fail("Deposited RMA scale inconsistent with metadata; investigate without outcome-dependent preprocessing")
registry<-fread("runs/001A_metadata/sample_registry_v1.tsv")
m<-registry[series=="GSE73661"];m<-m[match(colnames(x),gsm)]
stopifnot(identical(m$gsm,colnames(x)))
p<-as.data.frame(fread(plan$pair_file))
stopifnot(nrow(p)==23,sum(p$endoscopic_healing=="Yes")==8,!anyDuplicated(p$subject),all(p$baseline_gsm%in%colnames(x)),all(p$followup_gsm%in%colnames(x)))
p$endoscopic_healing<-factor(p$endoscopic_healing,levels=c("No","Yes"))
map<-fread("runs/001D_annotation/GPL6244_canonical_probe_mapping.tsv",na.strings="")
map<-map[match(rownames(x),probe)]
keep<-!is.na(map$retained)&map$retained==1&!is.na(map$canonical_symbol)
g<-avereps(x[keep,,drop=FALSE],ID=map$canonical_symbol[keep]);g<-g[order(rownames(g)),,drop=FALSE]
cand<-lock$candidate$id;infl<-"HALLMARK_INFLAMMATORY_RESPONSE"
cgenes<-lock$candidate$genes;common<-lock$inflammation_covariate$genes
if(!all(cgenes%in%rownames(g)))fail("Candidate six incomplete")
if(!all(common%in%rownames(g)))fail("Fixed common194 inflammation incomplete")
mem<-fread("runs/001D_annotation/program_membership_canonical.tsv")
refs<-lock$secondary_reference_programs;sets<-list();cov<-list()
for(s in c(cand,infl,refs)){
 genes<-unique(mem[program==s]$gene);hits<-intersect(genes,rownames(g));sets[[s]]<-hits
 cov[[length(cov)+1]]<-data.table(program=s,source_n=length(genes),measured_n=length(hits),coverage=length(hits)/length(genes),missing=paste(setdiff(genes,hits),collapse=";"))
}
stopifnot(setequal(sets[[cand]],cgenes))
sets[[infl]]<-common
S<-t(vapply(sets,function(z)colMeans(g[z,,drop=FALSE]),numeric(ncol(g))));colnames(S)<-colnames(g)
BL<-S[,p$baseline_gsm];PO<-S[,p$followup_gsm];D<-PO-BL;colnames(BL)<-colnames(PO)<-colnames(D)<-p$subject
saveRDS(list(gene_expression=g,probe_expression_log2=x,metadata=m,pairs=p,sample_scores=S,baseline=BL,post=PO,delta=D,sets=sets),file.path(out,"support_scores.rds"))
fwrite(rbindlist(cov),file.path(out,"actual_program_coverage.tsv"),sep="\t")
fwrite(data.table(gene=common,measured=common%in%rownames(g)),file.path(out,"scored_inflammation194.tsv"),sep="\t")
dd<-p;dd$delta_candidate<-D[cand,];dd$baseline_candidate<-BL[cand,];dd$delta_common194_inflammation<-D[infl,]
fwrite(dd,file.path(out,"support_patient_data.tsv"),sep="\t")
term<-"endoscopic_healingYes"
rowfit<-function(f,s,label,method="OLS"){
 X<-model.matrix(f);a<-summary(f)$coef[term,];df<-df.residual(f);be<-unname(a[1]);se<-unname(a[2])
 if(method=="HC3"){bread<-solve(crossprod(X));e<-resid(f)/(1-hatvalues(f));vv<-bread%*%crossprod(X,X*as.numeric(e^2))%*%bread;se<-sqrt(vv[term,term])}
 data.table(program=s,model=label,method=method,estimate=be,se=se,lower=be-qt(.975,df)*se,upper=be+qt(.975,df)*se,p=2*pt(-abs(be/se),df),n=nobs(f),df=df,rank=f$rank,standardized_condition=kappa(cbind(1,scale(X[,-1,drop=FALSE])),exact=TRUE),max_VIF=max(diag(solve(cor(X[,-1,drop=FALSE])))),max_leverage=max(hatvalues(f)),evaluable=TRUE)
}
forms<-plan$models
effects<-list();fits<-list();issues<-list()
for(label in names(forms)){
 f<-lm(as.formula(forms[[label]]),dd)
 if(f$rank!=length(coef(f))||!term%in%rownames(summary(f)$coef)){issues[[length(issues)+1]]<-data.table(model=label,reason="Rank deficient or healing coefficient not estimable");next}
 fits[[label]]<-f
 for(method in c("OLS","HC3"))effects[[length(effects)+1]]<-rowfit(f,cand,label,method)
}
if(!"adjusted_support"%in%names(fits))fail("Fixed adjusted support model not estimable")
ee<-rbindlist(effects);fwrite(ee,file.path(out,"candidate_support_models.tsv"),sep="\t")
saveRDS(fits,file.path(out,"support_model_fits.rds"))
f<-fits$adjusted_support
fwrite(data.table(subject=dd$subject,endoscopic_healing=dd$endoscopic_healing,leverage=hatvalues(f),cooks_distance=cooks.distance(f),residual=resid(f),studentized_residual=rstudent(f)),file.path(out,"patient_diagnostics.tsv"),sep="\t")
X<-model.matrix(f)[,-1,drop=FALSE];fwrite(data.table(term=colnames(X),VIF=diag(solve(cor(X)))),file.path(out,"adjusted_support_VIF.tsv"),sep="\t")
loo<-rbindlist(lapply(seq_len(nrow(dd)),function(i){
 ff<-lm(formula(f),dd[-i,]);ok<-ff$rank==length(coef(ff))&&term%in%rownames(summary(ff)$coef)
 data.table(omitted_subject=dd$subject[i],evaluable=ok,estimate=if(ok)unname(coef(ff)[term]) else NA_real_,p=if(ok)summary(ff)$coef[term,4] else NA_real_)
}))
fwrite(loo,file.path(out,"leave_one_patient_out.tsv"),sep="\t")
sec<-list()
for(s in refs){
 covs<-rbindlist(cov)[program==s]
 if(covs$coverage<.8){sec[[length(sec)+1]]<-data.table(program=s,evaluable=FALSE,p=NA_real_,reason="Coverage below80%");next}
 rr<-dd;rr$delta_candidate<-D[s,];rr$baseline_candidate<-BL[s,];ff<-lm(as.formula(forms$adjusted_support),rr)
 if(ff$rank!=length(coef(ff))){sec[[length(sec)+1]]<-data.table(program=s,evaluable=FALSE,p=NA_real_,reason="Rank deficient");next}
 sec[[length(sec)+1]]<-rowfit(ff,s,"reference_descriptive","OLS")
}
sr<-rbindlist(sec,fill=TRUE);sr$BH_across_three_descriptive<-p.adjust(sr$p,"BH",n=3);fwrite(sr,file.path(out,"reference_descriptive_models.tsv"),sep="\t")
# All cohorts remain separate. No pooled P value is calculated.
scoretab<-data.table(subject=dd$subject,endoscopic_healing=dd$endoscopic_healing,baseline_candidate=BL[cand,],post_candidate=PO[cand,],delta_candidate=D[cand,],delta_inflammation=D[infl,],baseline_endoscopic=dd$baseline_mayo_endoscopic,post_endoscopic=dd$followup_mayo_endoscopic)
fwrite(scoretab,file.path(out,"paired_score_plot_source.tsv"),sep="\t")
if(length(issues))fwrite(rbindlist(issues),file.path(out,"non_evaluable_auxiliary_models.tsv"),sep="\t")
d1<-rbindlist(cov);dups<-colnames(x)[duplicated(as.data.frame(t(x)))]
summary<-list(status="COMPLETE",evaluable=TRUE,role=plan$role,can_replace_GSE23597=FALSE,endpoint="Post-treatment Mayo endoscopic subscore0or1; W4/6 source window",n_pairs=nrow(p),healing=8,nonhealing=15,arrays=ncol(x),probes=nrow(x),canonical_genes=nrow(g),source_quantiles=as.list(q),normalization="Deposited aroma.affymetrix RMA log2 retained; no second normalization or log2",exact_duplicate_arrays=dups,patient_exclusions_by_expression=FALSE,adjusted_support=as.list(ee[model=="adjusted_support"&method=="OLS"]),residual_sigma=summary(f)$sigma,candidate_delta_SD=sd(D[cand,]),GSE23597_primary_support_remains_FALSE=TRUE,IBrD_not_compared=TRUE,pooled_P_computed=FALSE,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
write_json(summary,file.path(out,"support_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=15)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"));print(ee);print(sr);cat("GSE73661_SUPPORT_COMPLETE\n")
