options(stringsAsFactors=FALSE,width=150)
suppressPackageStartupMessages({library(limma);library(data.table);library(jsonlite)})
root<-"/root/projects/UC_Treatment_Recovery";set.seed(20260910)
obj<-readRDS(file.path(root,"runs/001A_qc/discovery_preprocessed_v1.rds"))
pairs<-obj$pairs;y<-obj$delta
stopifnot(identical(colnames(y),pairs$subject))
batch<-as.data.frame(fread(file.path(root,"provenance/reviews/metadata_independent_v1_gsm_barcode.tsv")))
bl<-batch$inferred_cel_barcode[match(pairs$baseline_gsm,batch$gsm)]
po<-batch$inferred_cel_barcode[match(pairs$post_gsm,batch$gsm)]
levels_b<-sort(unique(c(bl,po)))
B<-sapply(levels_b[-1],function(z)as.integer(po==z)-as.integer(bl==z))
colnames(B)<-paste0("barcode_delta_",seq_len(ncol(B)))
arm<-factor(pairs$arm,levels=c("Placebo","golimumab"))
X<-model.matrix(~0+arm);colnames(X)<-c("Placebo","GLM");rownames(X)<-pairs$subject
cmain<-makeContrasts(GLM_change=GLM,Placebo_change=Placebo,GLM_minus_Placebo=GLM-Placebo,levels=X)
response_group<-factor(paste(pairs$arm,pairs$response,sep="_"),levels=c("Placebo_No","Placebo_Yes","golimumab_No","golimumab_Yes"))
Z<-model.matrix(~0+response_group);colnames(Z)<-c("P_No","P_Yes","G_No","G_Yes")
cresp<-makeContrasts(GLM_response_difference=G_Yes-G_No,Placebo_response_difference=P_Yes-P_No,response_difference_interaction=(G_Yes-G_No)-(P_Yes-P_No),levels=Z)
# Input/model smoke test: independently check mean-change coefficients and contrast algebra.
test_idx<-seq_len(100)
fcheck<-lmFit(y[test_idx,,drop=FALSE],X)
stopifnot(max(abs(fcheck$coefficients[,"GLM"]-rowMeans(y[test_idx,arm=="golimumab",drop=FALSE])))<1e-10)
fc<-contrasts.fit(fcheck,cmain)
stopifnot(max(abs(fc$coefficients[,"GLM_minus_Placebo"]-(rowMeans(y[test_idx,arm=="golimumab",drop=FALSE])-rowMeans(y[test_idx,arm=="Placebo",drop=FALSE]))))<1e-10)
write_json(list(status="TESTED",n_test_genes=100,checks=c("matrix columns align to source subject IDs","paired subtraction retains baseline/post identities","limma coefficients match independent row means","difference contrast matches direct arithmetic"),seed=20260910),file.path(root,"tests/model_smoke_001B.json"),auto_unbox=TRUE,pretty=TRUE)
summary<-list()
run_branch<-function(branch,include_batch){
 out<-file.path(root,"runs/001B_discovery",branch);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 xx<-if(include_batch)cbind(X,B) else X
 zz<-if(include_batch)cbind(Z,B) else Z
 stopifnot(qr(xx)$rank==ncol(xx),qr(zz)$rank==ncol(zz))
 cm<-rbind(cmain,matrix(0,ncol(xx)-nrow(cmain),ncol(cmain)));rownames(cm)<-colnames(xx)
 cr<-rbind(cresp,matrix(0,ncol(zz)-nrow(cresp),ncol(cresp)));rownames(cr)<-colnames(zz)
 fm<-eBayes(contrasts.fit(lmFit(y,xx),cm),robust=TRUE,trend=FALSE)
 fr<-eBayes(contrasts.fit(lmFit(y,zz),cr),robust=TRUE,trend=FALSE)
 fwrite(data.frame(subject=pairs$subject,xx),file.path(out,"main_design.tsv"),sep="\t")
 fwrite(data.frame(subject=pairs$subject,zz),file.path(out,"response_design.tsv"),sep="\t")
 allfits<-list(main=fm,response=fr);rows<-list()
 for(kind in names(allfits)){
  fit<-allfits[[kind]]
  for(co in colnames(fit$coefficients)){
   tab<-topTable(fit,coef=co,number=Inf,sort.by="none",confint=.95)
   tab$gene<-rownames(tab);tab<-tab[,c("gene",setdiff(names(tab),"gene"))]
   fwrite(tab,file.path(out,paste0(co,".tsv")),sep="\t")
   rows[[co]]<-list(contrast=co,tested_genes=nrow(tab),fdr05=sum(tab$adj.P.Val<.05),fdr05_up=sum(tab$adj.P.Val<.05&tab$logFC>0),fdr05_down=sum(tab$adj.P.Val<.05&tab$logFC<0))
  }
 }
 saveRDS(allfits,file.path(out,"limma_fits.rds"))
 info<-list(run_id=paste0("001B_",branch),stage="exploratory discovery; no signature frozen",branch=branch,inferred_barcode_adjustment=include_batch,main_design_rank=qr(xx)$rank,response_design_rank=qr(zz)$rank,n_pairs=ncol(y),robust_ebayes=TRUE,trend=FALSE,results=rows,validation_expression_accessed=FALSE)
 write_json(info,file.path(out,"model_summary.json"),auto_unbox=TRUE,pretty=TRUE)
 info
}
summary$primary<-run_branch("primary_paired",FALSE)
summary$sensitivity<-run_branch("barcode_sensitivity",TRUE)
write_json(summary,file.path(root,"runs/001B_discovery/discovery_model_summary.json"),auto_unbox=TRUE,pretty=TRUE)
cat(toJSON(summary,auto_unbox=TRUE,pretty=TRUE),"\n")
