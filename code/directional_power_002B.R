options(stringsAsFactors=FALSE);suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery");out<-"runs/002B_precision"
a<-fread(file.path(out,"fitted_precision_components.tsv"));v<-a[cohort=="GSE23597"];d<-a[cohort=="GSE92415"]
pplus<-function(e)pt(qt(.975,v$df),v$df,ncp=e/v$SE,lower.tail=FALSE)
ptwo<-function(e)pplus(e)+pt(-qt(.975,v$df),v$df,ncp=e/v$SE)
eff<-c(0,.25,.4,d$estimate,.8,1,v$estimate)
tab<-data.table(assumed_true_effect=eff,two_sided_rejection_power=vapply(eff,ptwo,numeric(1)),positive_direction_support_probability=vapply(eff,pplus,numeric(1)))
tab$negative_direction_rejection_probability<-tab$two_sided_rejection_power-tab$positive_direction_support_probability
fwrite(tab,file.path(out,"directional_support_power_audit.tsv"),sep="\t")
summary<-list(two_sided_MDE80=uniroot(function(e)ptwo(e)-.8,c(0,10),tol=1e-11)$root,direction_consistent_support_MDE80=uniroot(function(e)pplus(e)-.8,c(0,10),tol=1e-11)$root,interpretation="Original tables give two-sided rejection probability. Frozen support also requires positive coefficient, whose conditional probability excludes negative-tail rejection. Neither is an observed true power or explanation of negative result.")
write_json(summary,file.path(out,"directional_power_definition.json"),auto_unbox=TRUE,pretty=TRUE,digits=15);print(tab)
