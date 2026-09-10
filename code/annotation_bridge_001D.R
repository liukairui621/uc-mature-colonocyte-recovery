options(stringsAsFactors=FALSE);suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd("/root/projects/UC_Treatment_Recovery")
out<-"runs/001D_annotation"
o<-readRDS(file.path(out,"discovery_preprocessed_canonical.rds"));s<-readRDS(file.path(out,"program_scores.rds"))
mem<-fread(file.path(out,"program_membership_canonical.tsv"));infl<-"HALLMARK_INFLAMMATORY_RESPONSE";cand<-"MATURE_COLONOCYTE6_EXPLORATORY"
genes<-mem[program==infl]$gene
full<-colMeans(o$gene_expression[genes,,drop=FALSE])
p<-s$pairs;p$delta<-s$delta[cand,];p$baseline<-s$baseline[cand,];p$inflam<-full[p$post_gsm]-full[p$baseline_gsm]
f<-lm(delta~response+arm+baseline+inflam,p);a<-summary(f)$coef["responseYes",];ci<-confint(f,"responseYes")
r<-data.table(annotation="HGNC canonical",inflammation_panel="Full 200 measurable in discovery",estimate=unname(a[1]),se=unname(a[2]),lower=ci[1],upper=ci[2],p=unname(a[4]))
fwrite(r,file.path(out,"annotation_only_full200_discovery_bridge.tsv"),sep="\t")
# Document exact algebraic equality of change and post-score models.
p$post_score<-p$delta+p$baseline
fp<-lm(post_score~response+arm+baseline+inflam,p)
stopifnot(abs(coef(fp)["responseYes"]-coef(f)["responseYes"])<1e-10)
write_json(list(delta_post_response_coefficient_max_error=abs(coef(fp)["responseYes"]-coef(f)["responseYes"]),meaning="Algebraic equivalence for identical covariates; does not confer causality"),file.path(out,"baseline_model_identity_check.json"),auto_unbox=TRUE,pretty=TRUE,digits=15)
print(r)
