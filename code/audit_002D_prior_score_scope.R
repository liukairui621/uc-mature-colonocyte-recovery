setwd("/root/projects/UC_Treatment_Recovery")
z<-readRDS("runs/002B_GSE73661_support/support_scores.rds")
m<-z$metadata
vdz<-grepl("^vdz",m$treatment);ifx<-m$treatment=="IFX"
stopifnot(all(m$gsm[vdz] %in% colnames(z$sample_scores)),nrow(z$pairs)==23,all(z$pairs$treatment_path=="IFX"))
a<-list(existing_object_names=names(z),sample_score_dimensions=dim(z$sample_scores),scored_VDZ_arrays=sum(vdz),scored_IFX_arrays=sum(ifx),paired_model_patients=nrow(z$pairs),paired_model_treatments=unique(z$pairs$treatment_path),prior_VDZ_sample_scores_computed=TRUE,values_printed=FALSE,prior_VDZ_association_model_found_in_retained_code=FALSE,caveat="No evidence found in retained project code; this is not proof of absence of all historical or external inspection")
jsonlite::write_json(a,"runs/002D_VDZ_metadata/prior_score_scope.json",pretty=TRUE,auto_unbox=TRUE)
print(a)
