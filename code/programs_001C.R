options(stringsAsFactors=FALSE,width=160);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-"/root/projects/UC_Treatment_Recovery";setwd(root)
out<-"runs/001C_programs";dir.create(out,showWarnings=FALSE,recursive=TRUE)
o<-readRDS("runs/001A_qc/discovery_preprocessed_v1.rds");g<-o$gene_expression;p<-as.data.frame(o$pairs)
stopifnot(identical(colnames(o$delta),p$subject))
p$arm<-factor(p$arm,levels=c("Placebo","golimumab"));p$response<-factor(p$response,levels=c("No","Yes"))
p$mayo_baseline<-as.numeric(p$mayo_baseline);p$mayo_post<-as.numeric(p$mayo_post)
mem<-fread("inputs/annotation/program_membership_001C.tsv",sep="\t")
stopifnot(!anyDuplicated(mem[,.(program,gene)]),all(is.finite(mem$weight)))
mem[,measured:=gene%in%rownames(g)]
coverage<-mem[,.(source_genes=.N,measured_genes=sum(measured),coverage=mean(measured),missing_genes=paste(gene[!measured],collapse=";")),by=program]
fwrite(coverage,file.path(out,"program_coverage.tsv"),sep="\t")
sets<-split(mem$gene,mem$program)
S<-t(vapply(names(sets),function(s){
 mm<-mem[program==s&measured==TRUE];stopifnot(nrow(mm)>=3)
 colSums(g[mm$gene,,drop=FALSE]*mm$weight)/sum(abs(mm$weight))
},numeric(ncol(g))))
rownames(S)<-names(sets);colnames(S)<-colnames(g)
BL<-S[,match(p$baseline_gsm,colnames(S)),drop=FALSE];PO<-S[,match(p$post_gsm,colnames(S)),drop=FALSE];D<-PO-BL
colnames(BL)<-colnames(PO)<-colnames(D)<-p$subject
fwrite(data.table(gsm=colnames(S),t(S)),file.path(out,"program_all_sample_scores.tsv"),sep="\t")
fwrite(data.table(subject=p$subject,arm=p$arm,response=p$response,t(D)),file.path(out,"program_patient_delta.tsv"),sep="\t")
fwrite(data.table(subject=p$subject,t(BL)),file.path(out,"program_patient_baseline.tsv"),sep="\t")
saveRDS(list(sample_scores=S,baseline=BL,post=PO,delta=D,pairs=p,membership=mem),file.path(out,"program_scores.rds"))
fam<-function(s)ifelse(grepl("^MARS_TableS1",s),"MARS_TableS1",ifelse(grepl("^MARS_MSigDB74",s),"MARS_MSigDB74","focused_and_reference_programs"))
ci<-function(fit,term){
 a<-summary(fit)$coefficients[term,];z<-confint(fit,term,level=.95)
 list(estimate=unname(a[1]),se=unname(a[2]),lower=unname(z[1]),upper=unname(z[2]),p=unname(a[4]),df=fit$df.residual,n=nobs(fit))
}
rowfit<-function(fit,term,program,model){
 as.data.table(c(list(program=program,model=model,term=term,family=fam(program)),ci(fit,term)))
}
scoredata<-function(s){
 dd<-p;dd$delta<-as.numeric(D[s,]);dd$baseline<-as.numeric(BL[s,])
 dd$inflam_delta<-as.numeric(D["HALLMARK_INFLAMMATORY_RESPONSE",]);dd$epi_delta<-as.numeric(D["EPITHELIAL_ABUNDANCE_PROXY4",])
 dd
}
batch<-fread("provenance/reviews/metadata_independent_v1_gsm_barcode.tsv",sep="\t")
bb<-batch$inferred_cel_barcode[match(p$baseline_gsm,batch$gsm)];bp<-batch$inferred_cel_barcode[match(p$post_gsm,batch$gsm)]
lv<-sort(unique(c(bb,bp)));B<-sapply(lv[-1],function(z)as.integer(bp==z)-as.integer(bb==z));colnames(B)<-paste0("barcode_delta_",seq_len(ncol(B)))
focus<-c("MATURE_COLONOCYTE6_EXPLORATORY","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
effects<-list();within<-list();fits<-list()
for(s in rownames(D)){
 dd<-scoredata(s)
 f<-lm(delta~arm+response,data=dd);effects[[length(effects)+1]]<-rowfit(f,"responseYes",s,"response_arm_adjusted");fits[[paste(s,"primary",sep="::")]]<-f
 f<-lm(delta~arm+response+baseline+mayo_baseline,data=dd);effects[[length(effects)+1]]<-rowfit(f,"responseYes",s,"response_baseline_adjusted")
 ddB<-cbind(dd,B);f<-lm(reformulate(c("arm","response",colnames(B)),"delta"),data=ddB)
 effects[[length(effects)+1]]<-rowfit(f,"responseYes",s,"response_barcode_adjusted")
 f<-lm(delta~arm*response,data=dd);effects[[length(effects)+1]]<-rowfit(f,"armgolimumab:responseYes",s,"arm_response_interaction")
 for(a in levels(p$arm)){
  sub<-dd[dd$arm==a,];f<-lm(delta~response,data=sub)
  effects[[length(effects)+1]]<-rowfit(f,"responseYes",s,paste0(a,"_response_difference"))
  z<-t.test(sub$delta)
  within[[length(within)+1]]<-data.table(program=s,arm=a,n=nrow(sub),mean_change=mean(sub$delta),lower=z$conf.int[1],upper=z$conf.int[2],p=z$p.value,unit="mean log2 post-minus-baseline")
 }
 f<-lm(delta~arm+baseline+mayo_baseline,data=dd);effects[[length(effects)+1]]<-rowfit(f,"armgolimumab",s,"arm_baseline_adjusted")
 if(s%in%focus){
  for(model in c("response_plus_inflammation","response_plus_epithelial_proxy","response_plus_both")){
   extra<-switch(model,response_plus_inflammation="inflam_delta",response_plus_epithelial_proxy="epi_delta",response_plus_both=c("inflam_delta","epi_delta"))
   f<-lm(reformulate(c("arm","response","baseline","mayo_baseline",extra),"delta"),data=dd)
   effects[[length(effects)+1]]<-rowfit(f,"responseYes",s,model);fits[[paste(s,model,sep="::")]]<-f
  }
 }
}
ee<-rbindlist(effects);stopifnot(all(is.finite(ee$estimate)),all(is.finite(ee$p)),all(ee$lower<=ee$upper))
ee[,FDR:=p.adjust(p,"BH"),by=.(model,family)]
fwrite(ee,file.path(out,"program_response_and_arm_effects_CI.tsv"),sep="\t")
ww<-rbindlist(within);ww[,family:=fam(program)];ww[,FDR:=p.adjust(p,"BH"),by=.(arm,family)]
fwrite(ww,file.path(out,"program_within_arm_change_CI.tsv"),sep="\t")
saveRDS(fits,file.path(out,"patient_model_fits.rds"))
# Compare candidate and OXPHOS with all prespecified reference sets, at the patient unit.
corrs<-list();core<-c("MATURE_COLONOCYTE6_EXPLORATORY","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
for(s in core)for(tar in setdiff(rownames(D),s)){
 ct<-cor.test(D[s,],D[tar,],method="pearson");sp<-cor.test(D[s,],D[tar,],method="spearman",exact=FALSE)
 ov<-intersect(sets[[s]],sets[[tar]])
 corrs[[length(corrs)+1]]<-data.table(program=s,reference=tar,n=ncol(D),pearson=unname(ct$estimate),lower=ct$conf.int[1],upper=ct$conf.int[2],p=ct$p.value,spearman=unname(sp$estimate),spearman_p=sp$p.value,shared_genes=length(ov),jaccard=length(ov)/length(union(sets[[s]],sets[[tar]])),shared_gene_symbols=paste(ov,collapse=";"))
}
cr<-rbindlist(corrs);cr[,FDR:=p.adjust(p,"BH"),by=program]
fwrite(cr,file.path(out,"program_overlap_and_delta_correlations.tsv"),sep="\t")
# Incremental contemporaneous association, not baseline prediction. Report estimates and intervals even when LRT is nonsignificant.
logrows<-list()
for(s in focus){
 dd<-scoredata(s);dd$resp_binary<-as.integer(dd$response=="Yes")
 for(spec in c("arm_inflammation","arm_inflammation_baseline")){
  vars<-c("arm","inflam_delta");if(spec=="arm_inflammation_baseline")vars<-c(vars,"baseline","mayo_baseline")
  f0<-glm(reformulate(vars,"resp_binary"),family=binomial(),data=dd)
  f1<-glm(reformulate(c(vars,"delta"),"resp_binary"),family=binomial(),data=dd)
  a<-summary(f1)$coef["delta",];lrt<-anova(f0,f1,test="LRT")
  be<-unname(a[1]);se<-unname(a[2]);sdx<-sd(dd$delta)
  logrows[[length(logrows)+1]]<-data.table(program=s,model=spec,n=nrow(dd),responders=sum(dd$resp_binary),beta_per_log2=be,SE=se,lower_beta=be-qnorm(.975)*se,upper_beta=be+qnorm(.975)*se,OR_per_SD=exp(be*sdx),OR_SD_lower=exp((be-qnorm(.975)*se)*sdx),OR_SD_upper=exp((be+qnorm(.975)*se)*sdx),wald_p=unname(a[4]),LRT_p=lrt$"Pr(>Chi)"[2],converged=f1$converged)
 }
}
lr<-rbindlist(logrows);lr[,FDR:=p.adjust(LRT_p,"BH"),by=model];fwrite(lr,file.path(out,"incremental_response_association.tsv"),sep="\t")
# No gene is removed after these influence analyses.
influence<-list();logo<-list()
for(s in focus){
 dd<-scoredata(s)
 for(i in seq_len(nrow(dd))){
  ff<-lm(delta~arm+response,data=dd[-i,])
  influence[[length(influence)+1]]<-data.table(program=s,omitted_subject=dd$subject[i],estimate=unname(coef(ff)["responseYes"]))
 }
}
s<-"MATURE_COLONOCYTE6_EXPLORATORY";mm<-mem[program==s&measured==TRUE]
gene_delta<-o$delta[mm$gene,,drop=FALSE];gg<-cor(t(gene_delta))
fwrite(data.table(gene=rownames(gg),gg),file.path(out,"candidate_gene_delta_correlations.tsv"),sep="\t")
for(z in mm$gene){
 dd<-scoredata(s);dd$delta<-colMeans(gene_delta[setdiff(mm$gene,z),,drop=FALSE])
 f<-lm(delta~arm+response,data=dd)
 logo[[length(logo)+1]]<-as.data.table(c(list(omitted_gene=z),ci(f,"responseYes")))
}
fwrite(rbindlist(influence),file.path(out,"leave_one_patient_out.tsv"),sep="\t")
fwrite(rbindlist(logo),file.path(out,"leave_one_candidate_gene_out.tsv"),sep="\t")
genefx<-list()
for(z in unique(c(mm$gene,"ABCG2","SLC16A1"))){
 dd<-p;dd$delta<-as.numeric(o$delta[z,]);f<-lm(delta~arm+response,data=dd)
 genefx[[length(genefx)+1]]<-as.data.table(c(list(gene=z,responder_mean_delta=mean(dd$delta[dd$response=="Yes"]),nonresponder_mean_delta=mean(dd$delta[dd$response=="No"])),ci(f,"responseYes")))
}
gf<-rbindlist(genefx);gf[,FDR_within_eight_selected:=p.adjust(p,"BH")]
genome<-fread("runs/001B2_calibration/primary_paired/pooled_response_genes.tsv",sep="\t")
gf[,genomewide_limma_FDR:=genome$adj.P.Val[match(gene,genome$gene)]]
fwrite(gf,file.path(out,"selected_gene_effects_exploratory.tsv"),sep="\t")
info<-list(status="COMPLETE",run="001C_programs",paired_patients=nrow(p),responders=sum(p$response=="Yes"),nonresponders=sum(p$response=="No"),programs=nrow(D),candidate="MATURE_COLONOCYTE6_EXPLORATORY",candidate_genes=mm$gene,candidate_source="Post hoc reviewer-suggested set; not claimed new genes or independently validated signature",score="unweighted mean log2 expression; pairwise delta",original_MARS_or_IBrD_score_reproduced=FALSE,validation_expression_accessed=FALSE,candidate_frozen=FALSE,notes=c("All regressions observational contemporaneous associations","MARS original-table and version-7.4 gene memberships are distinct prespecified sensitivity comparators","IBrD is a measured-gene mean proxy, not original PCA/pseudotime prediction","Epithelial4 is a coarse marker-expression proxy, not calibrated cell proportions","No negative association test interpreted as no effect or proof of low power"))
write_json(info,file.path(out,"program_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=10)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat(toJSON(info,auto_unbox=TRUE,pretty=TRUE),"\nPROGRAMS_COMPLETE\n")
