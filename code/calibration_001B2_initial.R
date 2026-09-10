options(stringsAsFactors=FALSE,width=160)
suppressPackageStartupMessages({library(limma);library(data.table);library(jsonlite)})
set.seed(20260910)
root<-"/root/projects/UC_Treatment_Recovery";setwd(root)
out<-file.path(root,"runs/001B2_calibration");dir.create(out,showWarnings=FALSE,recursive=TRUE)
o<-readRDS("runs/001A_qc/discovery_preprocessed_v1.rds");y<-o$delta;genes<-rownames(y);p<-o$pairs
stopifnot(identical(colnames(y),p$subject),ncol(y)==65)
p$arm<-factor(p$arm,levels=c("Placebo","golimumab"));p$response<-factor(p$response,levels=c("No","Yes"))
batch<-fread("provenance/reviews/metadata_independent_v1_gsm_barcode.tsv",sep="\t")
bl<-batch$inferred_cel_barcode[match(p$baseline_gsm,batch$gsm)]
po<-batch$inferred_cel_barcode[match(p$post_gsm,batch$gsm)]
lv<-sort(unique(c(bl,po)));B<-sapply(lv[-1],function(z)as.integer(po==z)-as.integer(bl==z))
colnames(B)<-paste0("barcode_delta_",seq_len(ncol(B)))
sets<-readRDS("inputs/annotation/gene_sets_001B.rds")$Hallmark
idx<-lapply(sets,function(g)which(genes%in%g));idx<-idx[lengths(idx)>=15&lengths(idx)<=500]
H<-t(vapply(idx,function(ii)colMeans(y[ii,,drop=FALSE]),numeric(ncol(y))))
rownames(H)<-names(idx);colnames(H)<-p$subject
fwrite(data.table(subject=p$subject,arm=p$arm,response=p$response,t(H)),file.path(out,"Hallmark_patient_delta.tsv"),sep="\t")
ci_t<-function(x){
 n<-length(x);est<-mean(x);se<-sd(x)/sqrt(n);c(estimate=est,se=se,lower=est-qt(.975,n-1)*se,upper=est+qt(.975,n-1)*se,p=2*pt(-abs(est/se),n-1),n=n)
}
ci_lm<-function(fit,term){
 cc<-summary(fit)$coefficients[term,];ci<-confint(fit,term,level=.95)
 c(estimate=cc[1],se=cc[2],lower=ci[1],upper=ci[2],p=cc[4],df=fit$df.residual)
}
eff<-list()
for(s in rownames(H)){
 dd<-as.data.frame(p);dd$delta<-H[s,];w<-t.test(dd$delta[dd$arm=="golimumab"],dd$delta[dd$arm=="Placebo"])
 for(a in levels(p$arm))eff[[length(eff)+1]]<-data.table(set=s,contrast=paste0(a,"_within_change"),as.list(ci_t(dd$delta[dd$arm==a])),method="one-sample t",unit="mean log2 post-minus-baseline")
 eff[[length(eff)+1]]<-data.table(set=s,contrast="GLM_minus_Placebo",estimate=unname(w$estimate[1]-w$estimate[2]),se=unname(w$stderr),lower=w$conf.int[1],upper=w$conf.int[2],p=w$p.value,df=unname(w$parameter),method="Welch t",unit="difference of mean log2 changes")
 for(a in levels(p$arm)){
  sub<-dd[dd$arm==a,];fit<-lm(delta~response,data=sub)
  eff[[length(eff)+1]]<-data.table(set=s,contrast=paste0(a,"_response_difference"),as.list(ci_lm(fit,"responseYes")),method="OLS t",unit="difference of mean log2 changes")
 }
 fit<-lm(delta~arm+response,data=dd)
 eff[[length(eff)+1]]<-data.table(set=s,contrast="pooled_response_arm_adjusted",as.list(ci_lm(fit,"responseYes")),method="OLS t",unit="arm-adjusted response difference in log2 change")
 fit<-lm(delta~arm*response,data=dd)
 eff[[length(eff)+1]]<-data.table(set=s,contrast="arm_response_interaction",as.list(ci_lm(fit,"armgolimumab:responseYes")),method="OLS t",unit="difference of response differences")
}
effects<-rbindlist(eff,fill=TRUE);effects[,FDR:=p.adjust(p,"BH"),by=contrast]
fwrite(effects,file.path(out,"Hallmark_patient_effects_CI.tsv"),sep="\t")
all<-list();gsumm<-list()
for(branch in c("primary_paired","barcode_sensitivity")){
 bdir<-file.path(out,branch);dir.create(bdir,showWarnings=FALSE)
 X0<-as.matrix(fread(file.path("runs/001B_discovery",branch,"main_design.tsv"),sep="\t")[,-1])
 Z0<-as.matrix(fread(file.path("runs/001B_discovery",branch,"response_design.tsv"),sep="\t")[,-1])
 W0<-model.matrix(~arm+response,data=p)
 if(branch=="barcode_sensitivity")W0<-cbind(W0,B)
 stopifnot(qr(W0)$rank==ncol(W0))
 fwrite(data.table(subject=p$subject,W0),file.path(bdir,"pooled_design.tsv"),sep="\t")
 cc<-rep(0,ncol(W0));cc[match("responseYes",colnames(W0))]<-1
 fit<-eBayes(contrasts.fit(lmFit(y,W0),cc),robust=TRUE,trend=FALSE)
 tt<-topTable(fit,coef=1,number=Inf,sort.by="none",confint=.95);tt$gene<-rownames(tt)
 names(tt)[names(tt)=="AveExpr"]<-"mean_patient_delta"
 names(tt)[names(tt)=="logFC"]<-"response_contrast_log2"
 fwrite(tt,file.path(bdir,"pooled_response_genes.tsv"),sep="\t");saveRDS(fit,file.path(bdir,"pooled_response_fit.rds"))
 gsumm[[branch]]<-list(fdr05=sum(tt$adj.P.Val<.05),up=sum(tt$adj.P.Val<.05&tt$response_contrast_log2>0),down=sum(tt$adj.P.Val<.05&tt$response_contrast_log2<0),design_rank=qr(W0)$rank)
 for(co in c("GLM_change","Placebo_change","GLM_minus_Placebo","GLM_response_difference","pooled_response_arm_adjusted")){
  xx<-if(co=="pooled_response_arm_adjusted")W0 else if(co=="GLM_response_difference")Z0 else X0
  cc<-setNames(rep(0,ncol(xx)),colnames(xx))
  if(co=="GLM_change")cc["GLM"]<-1
  if(co=="Placebo_change")cc["Placebo"]<-1
  if(co=="GLM_minus_Placebo")cc[c("GLM","Placebo")]<-c(1,-1)
  if(co=="GLM_response_difference")cc[c("G_Yes","G_No")]<-c(1,-1)
  if(co=="pooled_response_arm_adjusted")cc["responseYes"]<-1
  for(met in c("camera_fixed001","camera_estimated","camera_ranks001","fry")){
   cat(branch,co,met,"\n");flush.console()
   if(met=="fry"){
    z<-fry(y,index=idx,design=xx,contrast=cc,sort="none",standardize="residual.sd")
    z<-data.table(set=rownames(z),z);z$null_hypothesis<-"self-contained; directional P plus mixed P"
   }else{
    z<-camera(y,index=idx,design=xx,contrast=cc,sort=FALSE,inter.gene.cor=if(met=="camera_estimated")NA else .01,use.ranks=met=="camera_ranks001")
    z<-data.table(set=rownames(z),z);z$null_hypothesis<-"competitive vs other assayed genes"
   }
   z$branch<-branch;z$contrast<-co;z$method<-met;z$rank<-rank(z$PValue,ties.method="min")
   md<-file.path(bdir,met);dir.create(md,showWarnings=FALSE)
   fwrite(z,file.path(md,paste0(co,".tsv")),sep="\t");all[[length(all)+1]]<-z
  }
 }
}
res<-rbindlist(all,fill=TRUE);res[,FDR_across_five_contrasts:=p.adjust(PValue,"BH"),by=.(branch,method)]
fwrite(res,file.path(out,"all_Hallmark_tests.tsv"),sep="\t")
summary<-res[,.(sets=.N,FDR05=sum(FDR<.05),median_estimated_correlation=if(all(is.na(Correlation)))NA_real_ else median(Correlation,na.rm=TRUE)),by=.(branch,contrast,method)]
fwrite(summary,file.path(out,"method_sensitivity_counts.tsv"),sep="\t")
ranks<-list()
for(b in unique(res$branch))for(co in unique(res$contrast)){
 zz<-res[branch==b&contrast==co]
 for(met in c("camera_estimated","camera_ranks001","fry")){
  a<-zz[method=="camera_fixed001"];d<-zz[method==met];d<-d[match(a$set,d$set)]
  ranks[[length(ranks)+1]]<-data.table(branch=b,contrast=co,comparison=paste("camera_fixed001",met,sep="_"),spearman_rank=cor(a$rank,d$rank,method="spearman"),direction_agreement=mean(a$Direction==d$Direction))
 }
}
fwrite(rbindlist(ranks),file.path(out,"rank_direction_stability.tsv"),sep="\t")
tG<-fread("runs/001B_discovery/primary_paired/GLM_change.tsv",sep="\t")
tP<-fread("runs/001B_discovery/primary_paired/Placebo_change.tsv",sep="\t");tP<-tP[match(tG$gene,tP$gene)]
sel<-tG$adj.P.Val<.05
tD<-fread("runs/001B_discovery/primary_paired/GLM_minus_Placebo.tsv",sep="\t")
audit<-list(genome_effect_correlation=cor(tG$logFC,tP$logFC),selected_GLM_genes_effect_correlation=cor(tG$logFC[sel],tP$logFC[sel]),selected_GLM_genes_same_sign=mean(sign(tG$logFC[sel])==sign(tP$logFC[sel])),between_arm_abs_effect_95_quantile=unname(quantile(abs(tD$logFC),.95)),between_arm_median_gene_CI_width=median(tD$CI.R-tD$CI.L),pooled_response_genes=gsumm)
write_json(audit,file.path(out,"review_numeric_audit.json"),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat(toJSON(audit,pretty=TRUE,auto_unbox=TRUE),"\nCALIBRATION_COMPLETE\n")
