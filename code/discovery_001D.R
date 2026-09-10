options(stringsAsFactors=FALSE,width=160);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(limma);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/001D_annotation"
o<-readRDS("runs/001A_qc/discovery_preprocessed_v1.rds")
map<-fread(file.path(out,"GPL13158_canonical_probe_mapping.tsv"),na.strings="")
map<-map[match(rownames(o$probe_expression),map$probe)]
stopifnot(identical(map$probe,rownames(o$probe_expression)))
keep<-map$retained==1&!is.na(map$canonical_symbol)
g<-avereps(o$probe_expression[keep,,drop=FALSE],ID=map$canonical_symbol[keep])
g<-g[order(rownames(g)),,drop=FALSE];p<-as.data.frame(o$pairs)
p$arm<-factor(p$arm,levels=c("Placebo","golimumab"));p$response<-factor(p$response,levels=c("No","Yes"))
p$mayo_baseline<-as.numeric(p$mayo_baseline)
delta<-g[,p$post_gsm]-g[,p$baseline_gsm];colnames(delta)<-p$subject
cand<-"MATURE_COLONOCYTE6_EXPLORATORY";infl<-"HALLMARK_INFLAMMATORY_RESPONSE";CT<-"UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy"
mem<-fread(file.path(out,"program_membership_canonical.tsv"))
common<-fread(file.path(out,"frozen_common_inflammation_members.tsv"))$gene
stopifnot(all(common%in%rownames(g)))
source_mem<-copy(mem)
mem<-mem[!(program==infl&!gene%in%common)]
fwrite(mem,file.path(out,"program_scoring_membership.tsv"),sep="\t")
mem$measured<-mem$gene%in%rownames(g)
sets<-split(mem$gene,mem$program)
S<-t(vapply(names(sets),function(s){mm<-mem[program==s&measured];colSums(g[mm$gene,,drop=FALSE]*mm$weight)/sum(abs(mm$weight))},numeric(ncol(g))))
rownames(S)<-names(sets);colnames(S)<-colnames(g)
BL<-S[,p$baseline_gsm];PO<-S[,p$post_gsm];D<-PO-BL;colnames(BL)<-colnames(PO)<-colnames(D)<-p$subject
old<-readRDS("runs/001C_programs/program_scores.rds")
candidate_diff<-max(abs(S[cand,]-old$sample_scores[cand,]))
stopifnot(candidate_diff<1e-10)
saveRDS(list(gene_expression=g,delta=delta,pairs=p,metadata=o$metadata,annotation=map),file.path(out,"discovery_preprocessed_canonical.rds"))
saveRDS(list(sample_scores=S,baseline=BL,post=PO,delta=D,pairs=p,membership=mem),file.path(out,"program_scores.rds"))
fwrite(data.table(subject=p$subject,arm=p$arm,response=p$response,t(BL)),file.path(out,"program_baseline.tsv"),sep="\t")
fwrite(data.table(subject=p$subject,arm=p$arm,response=p$response,t(D)),file.path(out,"program_delta.tsv"),sep="\t")
scoredata<-function(s){dd<-p;dd$delta<-D[s,];dd$baseline<-BL[s,];dd$inflam<-D[infl,];dd$epi<-D["EPITHELIAL_ABUNDANCE_PROXY4",];dd$ct_delta<-D[CT,];dd$ct_baseline<-BL[CT,];dd}
rowfit<-function(fit,s,model,term="responseYes"){
 X<-model.matrix(fit);Z<-cbind(1,scale(X[,-1,drop=FALSE]));a<-summary(fit)$coef[term,];ci<-confint(fit,term)
 VIF<-diag(solve(cor(X[,-1,drop=FALSE])))
 data.table(program=s,model=model,estimate=unname(a[1]),se=unname(a[2]),lower=ci[1],upper=ci[2],p=unname(a[4]),n=nobs(fit),df=fit$df.residual,rank=fit$rank,raw_condition=kappa(X,exact=TRUE),standardized_condition=kappa(Z,exact=TRUE),max_VIF=max(VIF),max_leverage=max(hatvalues(fit)))
}
hc3<-function(fit,s,model){
 X<-model.matrix(fit);h<-hatvalues(fit);e<-resid(fit);bread<-solve(crossprod(X));v<-bread%*%crossprod(X,X*as.numeric((e/(1-h))^2))%*%bread
 be<-coef(fit)["responseYes"];se<-sqrt(v["responseYes","responseYes"]);df<-df.residual(fit)
 data.table(program=s,model=model,estimate=unname(be),se=unname(se),lower=unname(be-qt(.975,df)*se),upper=unname(be+qt(.975,df)*se),p=unname(2*pt(-abs(be/se),df)),n=nobs(fit),df=df)
}
focus<-c(cand,CT,"UC_IA13_JAKSTAT_2024_gene_mean_proxy","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
batch<-fread("provenance/reviews/metadata_independent_v1_gsm_barcode.tsv")
bb<-batch$inferred_cel_barcode[match(p$baseline_gsm,batch$gsm)];bp<-batch$inferred_cel_barcode[match(p$post_gsm,batch$gsm)]
lv<-sort(unique(c(bb,bp)));B<-sapply(lv[-1],function(z)as.integer(bp==z)-as.integer(bb==z));colnames(B)<-paste0("barcode",seq_len(ncol(B)))
rows<-list();robust<-list();fits<-list();logrows<-list();lopo<-list();diagrows<-list()
for(s in focus){
 dd<-cbind(scoredata(s),B)
 models<-list(arm_response=c("arm","response"),no_baseline_inflammation=c("arm","response","inflam"),baseline_only=c("arm","response","baseline"),primary_baseline_inflammation=c("arm","response","baseline","inflam"),baseline_inflammation_Mayo=c("arm","response","baseline","inflam","mayo_baseline"),baseline_inflammation_epithelial=c("arm","response","baseline","inflam","epi"),baseline_inflammation_barcode=c("arm","response","baseline","inflam",colnames(B)))
 if(s==cand){
  for(src in c("MARS_TableS1_column_","MARS_MSigDB74_")){
   sn<-grep(paste0("^",src),rownames(D),value=TRUE);mn<-paste0("M",seq_along(sn))
   for(j in seq_along(mn))dd[[paste0(src,mn[j])]]<-D[sn[j],]
   models[[paste0("baseline_inflammation_",src)]]<-c("arm","response","baseline","inflam",paste0(src,mn))
  }
  models$CT_conditioned<-c("arm","response","baseline","inflam","ct_delta","ct_baseline")
 }
 for(mod in names(models)){
  fit<-lm(reformulate(models[[mod]],"delta"),data=dd)
  stopifnot(fit$rank==length(coef(fit)))
  rows[[length(rows)+1]]<-rowfit(fit,s,mod);fits[[paste(s,mod,sep="::")]]<-fit
  if(mod%in%c("primary_baseline_inflammation","no_baseline_inflammation")){
   robust[[length(robust)+1]]<-hc3(fit,s,mod)
   if(s==cand)for(i in seq_len(nrow(dd)))lopo[[length(lopo)+1]]<-data.table(model=mod,omitted_subject=dd$subject[i],estimate=unname(coef(lm(formula(fit),data=dd[-i,]))["responseYes"]))
  }
  if(mod=="CT_conditioned"){
   X<-model.matrix(fit)[,-1,drop=FALSE]
   diagrows[[1]]<-data.table(term=colnames(X),VIF=diag(solve(cor(X))))
  }
 }
 if(s%in%c(cand,"HALLMARK_OXIDATIVE_PHOSPHORYLATION")){
  dd$resp<-as.integer(dd$response=="Yes")
  for(mod in c("no_baseline","own_baseline","own_baseline_and_Mayo")){
   vars<-c("arm","inflam")
   if(mod!="no_baseline")vars<-c(vars,"baseline")
   if(mod=="own_baseline_and_Mayo")vars<-c(vars,"mayo_baseline")
   f0<-glm(reformulate(vars,"resp"),binomial(),data=dd);f1<-glm(reformulate(c(vars,"delta"),"resp"),binomial(),data=dd)
   a<-summary(f1)$coef["delta",];lrt<-anova(f0,f1,test="LRT")
   logrows[[length(logrows)+1]]<-data.table(program=s,model=mod,beta=unname(a[1]),se=unname(a[2]),lower=unname(a[1]-qnorm(.975)*a[2]),upper=unname(a[1]+qnorm(.975)*a[2]),LRT_p=lrt$"Pr(>Chi)"[2],converged=f1$converged)
  }
 }
}
ee<-rbindlist(rows);ee$exploratory_FDR_within_model<-NA_real_
for(mod in unique(ee$model)){ix<-which(ee$model==mod);ee$exploratory_FDR_within_model[ix]<-p.adjust(ee$p[ix],"BH")}
fwrite(ee,file.path(out,"discovery_response_models.tsv"),sep="\t")
fwrite(rbindlist(robust),file.path(out,"discovery_HC3.tsv"),sep="\t")
fwrite(rbindlist(diagrows),file.path(out,"CT_conditioned_VIF.tsv"),sep="\t")
fwrite(rbindlist(lopo),file.path(out,"leave_one_patient_out.tsv"),sep="\t")
ll<-rbindlist(logrows);ll$FDR<-NA_real_
for(mod in unique(ll$model)){ix<-which(ll$model==mod);ll$FDR[ix]<-p.adjust(ll$LRT_p[ix],"BH")}
fwrite(ll,file.path(out,"incremental_logistic_association.tsv"),sep="\t")
saveRDS(fits,file.path(out,"discovery_model_fits.rds"))
cr<-list()
for(s in setdiff(rownames(D),cand)){
 ct<-cor.test(D[cand,],D[s,]);ov<-intersect(sets[[cand]],sets[[s]])
 cr[[length(cr)+1]]<-data.table(reference=s,r=unname(ct$estimate),lower=ct$conf.int[1],upper=ct$conf.int[2],p=ct$p.value,shared_n=length(ov),shared_genes=paste(ov,collapse=";"))
}
fwrite(rbindlist(cr),file.path(out,"candidate_reference_correlations.tsv"),sep="\t")
logo<-list()
for(z in sets[[cand]]){
 dd<-scoredata(cand);genes<-setdiff(sets[[cand]],z)
 dd$delta<-colMeans(delta[genes,,drop=FALSE]);dd$baseline<-colMeans(g[genes,p$baseline_gsm,drop=FALSE])
 for(mod in c("arm_response","primary_baseline_inflammation")){
  f<-if(mod=="arm_response")lm(delta~arm+response,data=dd) else lm(delta~arm+response+baseline+inflam,data=dd)
  a<-rowfit(f,cand,mod);a$omitted_gene<-z;logo[[length(logo)+1]]<-a
 }
}
fwrite(rbindlist(logo),file.path(out,"leave_one_gene_out.tsv"),sep="\t")
# Updated arm changes and between-arm contrasts for all scored programs.
ar<-list()
for(s in rownames(D)){
 dd<-scoredata(s);f<-lm(delta~arm,data=dd)
 a<-rowfit(f,s,"GLM_minus_Placebo","armgolimumab");ar[[length(ar)+1]]<-a
 for(arm in levels(p$arm)){z<-t.test(dd$delta[dd$arm==arm]);ar[[length(ar)+1]]<-data.table(program=s,model=paste0(arm,"_within"),estimate=unname(z$estimate),lower=z$conf.int[1],upper=z$conf.int[2],p=z$p.value,n=sum(dd$arm==arm),df=unname(z$parameter))}
}
fwrite(rbindlist(ar,fill=TRUE),file.path(out,"program_arm_changes.tsv"),sep="\t")
write_json(list(status="COMPLETE",patients=nrow(p),genes=nrow(g),candidate_score_max_abs_change=candidate_diff,inflammation_common_genes=length(common),primary_discovery_effect=as.list(ee[program==cand&model=="primary_baseline_inflammation"]),validation_expression_accessed=FALSE,IBrD_excluded=TRUE,all_discovery_tests_exploratory=TRUE),file.path(out,"discovery_amendment_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=15)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(ee[program==cand|model=="primary_baseline_inflammation"]);cat("001D_DISCOVERY_COMPLETE\n")
