options(stringsAsFactors=FALSE,width=150);suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/002B_precision"
m<-as.data.frame(fread("runs/001A_metadata/sample_registry_v1.tsv",na.strings=""));m<-m[m$series=="GSE73661"&m$treatment=="IFX",]
print(m[,c("subject","gsm","time","week","response","mayo_endoscopic")]);print(names(m))
