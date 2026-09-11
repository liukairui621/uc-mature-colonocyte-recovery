suppressPackageStartupMessages(library(jsonlite))
root<-normalizePath(".")
out<-file.path(root,"runs/003E_ambient_negative_controls")
exactp<-function(y,g){
 n<-length(y); k<-sum(g==1); obs<-mean(y[g==1])-mean(y[g==0])
 cmb<-combn(n,k); dif<-apply(cmb,2,function(ix) mean(y[ix])-mean(y[-ix]))
 mean(abs(dif)>=abs(obs)-1e-12)
}
verify_models<-function(patient_path,model_path){
 d<-read.delim(patient_path,check.names=FALSE)
 m<-read.delim(model_path,check.names=FALSE)
 vals<-setdiff(names(d),c("patient","remission"))
 rows<-lapply(vals,function(v){
  z<-d[,c("remission",v)]; names(z)[2]<-"y"; z<-z[is.finite(z$y),]
  z$g<-as.integer(z$remission=="Remission")
  fit<-lm(y~g,z); ci<-confint(fit)["g",]
  r<-m[m$metric==v,]
  data.frame(metric=v,estimate_diff=unname(coef(fit)["g"]-r$estimate),
   lower_diff=unname(ci[1]-r$lower),upper_diff=unname(ci[2]-r$upper),
   exact_diff=exactp(z$y,z$g)-r$p_exact_permutation)
 })
 do.call(rbind,rows)
}
neg<-verify_models(file.path(out,"negative_control_CT_patient_deltas.tsv"),file.path(out,"negative_control_CT_models.tsv"))
con<-verify_models(file.path(out,"CT_contamination_patient_deltas.tsv"),file.path(out,"CT_contamination_patient_models.tsv"))
nm<-read.delim(file.path(out,"negative_control_CT_models.tsv"),check.names=FALSE)
ix<-nm$metric%in%c("delta_PTPRC_logCPM","delta_COL1A1_logCPM")
q<-p.adjust(nm$p_exact_permutation[ix],"BH",n=2)
qdiff<-max(abs(q-nm$BH_negative_genes_n2[ix]))
dc<-read.delim("runs/003D_decontx/sample_state_decontx_all.tsv",check.names=FALSE)
dc<-dc[dc$final_analysis=="Non ileal CT colonocyte",]
cross<-read.delim(file.path(out,"CT_contamination_cross_sectional_sample_state.tsv"),check.names=FALSE)
wd<-sapply(c("contamination_mean","contamination_median","contamination_q90"),function(v){
 a<-dc[dc$Remission_status=="Remission",v]; b<-dc[dc$Remission_status=="Non_Remission",v]
 p<-wilcox.test(a,b,exact=FALSE,correct=TRUE)$p.value
 p-cross$wilcoxon_rank_sum_p[cross$metric==v]
})
mx<-max(abs(c(unlist(neg[,-1]),unlist(con[,-1]),qdiff,wd)))
res<-list(status=if(mx<1e-12)"PASS" else "FAIL",
 independent_implementation="R lm/confint, exhaustive label combinations, p.adjust and wilcox.test",
 maximum_absolute_difference=mx,negative_control_model_differences=neg,
 contamination_patient_model_differences=con,BH_n2_max_difference=qdiff,
 cross_sectional_wilcoxon_differences=as.list(wd))
write_json(res,file.path(out,"arithmetic_verification_R.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat(toJSON(res,auto_unbox=TRUE,pretty=TRUE,digits=16))
if(mx>=1e-12) quit(status=1)
