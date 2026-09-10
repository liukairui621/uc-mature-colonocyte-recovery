options(width=150);setwd("/root/projects/UC_Treatment_Recovery")
a<-readRDS("runs/001D_annotation/discovery_model_fits.rds");b<-readRDS("runs/002A_GSE23597/candidate_model_fits.rds")
print(names(a));print(names(b));print(names(model.frame(a[["MATURE_COLONOCYTE6_EXPLORATORY::primary_baseline_inflammation"]])));print(names(model.frame(b[["primary_baseline_inflammation"]])))
