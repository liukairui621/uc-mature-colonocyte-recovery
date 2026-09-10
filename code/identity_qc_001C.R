options(stringsAsFactors=FALSE,width=150);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-"/root/projects/UC_Treatment_Recovery";setwd(root)
out<-"runs/001C_programs/identity_qc";dir.create(out,recursive=TRUE,showWarnings=FALSE)
o<-readRDS("runs/001A_qc/discovery_preprocessed_v1.rds")
m<-as.data.table(o$metadata);m[,center:=sub("-.*","",subject)]
g<-o$gene_expression;x<-o$probe_expression
stopifnot(identical(m$gsm,colnames(g)),identical(m$gsm,colnames(x)))
for(sym in c("XIST","RPS4Y1","KDM5D"))m[[sym]]<-as.numeric(g[sym,])
m[,Y_marker_mean:=(RPS4Y1+KDM5D)/2]
fwrite(m[,.(gsm,subject,center,treatment,time,XIST,RPS4Y1,KDM5D,Y_marker_mean)],file.path(out,"sex_marker_expression.tsv"),sep="\t")
cc<-cor(x,method="pearson")
base<-which(m$time=="Week 0"&m$disease=="Ulcerative Colitis (UC)")
post<-which(m$time=="Week 6"&m$disease=="Ulcerative Colitis (UC)")
comparisons<-list()
for(j in post){
 tab<-data.table(post_gsm=m$gsm[j],post_subject=m$subject[j],baseline_gsm=m$gsm[base],baseline_subject=m$subject[base],probe_pearson=cc[base,j],same_center=m$center[base]==m$center[j],same_arm=m$treatment[base]==m$treatment[j],source_ID_match=m$subject[base]==m$subject[j],abs_XIST_difference=abs(m$XIST[base]-m$XIST[j]),abs_Y_marker_difference=abs(m$Y_marker_mean[base]-m$Y_marker_mean[j]))
 tab[,correlation_rank_all:=rank(-probe_pearson,ties.method="min")]
 tab[,correlation_rank_same_arm:=NA_real_]
 tab[same_arm==TRUE,correlation_rank_same_arm:=rank(-probe_pearson,ties.method="min")]
 tab[,correlation_rank_same_center_arm:=NA_real_]
 tab[same_arm==TRUE&same_center==TRUE,correlation_rank_same_center_arm:=rank(-probe_pearson,ties.method="min")]
 comparisons[[length(comparisons)+1]]<-tab
}
all<-rbindlist(comparisons)
known<-all[source_ID_match==TRUE]
p<-o$pairs
unknown<-all[!post_gsm%in%p$post_gsm&!baseline_gsm%in%p$baseline_gsm]
fwrite(all,file.path(out,"all_pre_post_probe_correlations.tsv"),sep="\t")
fwrite(known,file.path(out,"known_pair_diagnostic.tsv"),sep="\t")
fwrite(unknown[order(post_gsm,-same_arm,-same_center,-probe_pearson)],file.path(out,"unmatched_candidate_diagnostic_only.tsv"),sep="\t")
report<-list(status="COMPLETE",n_confirmed_pairs=nrow(known),n_unmatched_post=length(unique(unknown$post_gsm)),n_unmatched_baseline=length(unique(unknown$baseline_gsm)),known_pair_top_rank_all=sum(known$correlation_rank_all==1),known_pair_top_rank_same_arm=sum(known$correlation_rank_same_arm==1),known_pair_top_rank_same_center_arm=sum(known$correlation_rank_same_center_arm==1),known_pair_cor_quantiles=as.list(quantile(known$probe_pearson,c(0,.25,.5,.75,1))),other_pair_cor_quantiles=as.list(quantile(all[source_ID_match==FALSE]$probe_pearson,c(0,.25,.5,.75,1))),additional_pairs_accepted=0,interpretation="Sex marker expression and bulk-expression correlation are QC diagnostics, not genotype identity. No pair was reconstructed. Known-pair ranks quantify how poorly/strongly nearest correlation recovers documented identity in this experiment.")
write_json(report,file.path(out,"identity_qc_summary.json"),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat(toJSON(report,auto_unbox=TRUE,pretty=TRUE),"\nIDENTITY_QC_COMPLETE\n")
