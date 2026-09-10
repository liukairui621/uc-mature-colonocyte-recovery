options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(jsonlite))
setwd("/root/projects/UC_Treatment_Recovery")
cand<-"MATURE_COLONOCYTE6_EXPLORATORY"
d<-readRDS("runs/001D_annotation/discovery_model_fits.rds")
v<-readRDS("runs/002A_GSE23597/candidate_model_fits.rds")
fitlist<-list(GSE92415=d[[paste(cand,"primary_baseline_inflammation",sep="::")]],GSE23597=v$primary_baseline_inflammation,discovery_before_inflammation=d[[paste(cand,"baseline_only",sep="::")]],validation_before_inflammation=v$without_inflammation)
models<-lapply(fitlist,function(f){X<-model.matrix(f);list(columns=colnames(X),X=unname(X),y=as.numeric(model.response(model.frame(f))),response_column="responseYes",formula=paste(deparse(formula(f)),collapse=" "))})
od<-readRDS("runs/001D_annotation/discovery_preprocessed_canonical.rds");ov<-readRDS("runs/002A_GSE23597/validation_scores.rds")
genes<-c("AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12")
gd<-od$delta[genes,,drop=FALSE];gv<-ov$gene_expression[genes,ov$pairs$week8_gsm,drop=FALSE]-ov$gene_expression[genes,ov$pairs$baseline_gsm,drop=FALSE]
write_json(list(models=models,genes=genes,gene_deltas=list(discovery=unname(gd),validation=unname(gv))),"provenance/reviews/002B_precision_recompute_inputs_v1.json",pretty=TRUE,auto_unbox=TRUE,digits=17)
cat("INDEPENDENT_SOURCE_FIT_EXTRACTION_COMPLETE\n")
