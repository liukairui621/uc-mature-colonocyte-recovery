args<-commandArgs(trailingOnly=TRUE)
agg_path<-if(length(args)>=1) args[1] else "runs/003D_decontx/sample_state_decontx_all.tsv"
pair_path<-if(length(args)>=2) args[2] else "runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv"
raw_path<-if(length(args)>=3) args[3] else "runs/003_cell_context/fixed_epithelial_state_models.tsv"
outdir<-if(length(args)>=4) args[4] else "runs/003D_decontx_state_family"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
genes<-c("AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12")
a<-read.delim(agg_path,check.names=FALSE,stringsAsFactors=FALSE)
pairs<-read.delim(pair_path,check.names=FALSE,stringsAsFactors=FALSE)
a<-a[a$major=="Non_ileal_epithelium",]
for(g in genes) a[[paste0(g,"_logCPM")]]<-log2((a[[paste0(g,"_decontx_counts")]]+0.5)/(a$decontx_total_UMI+1)*1e6)
a$state_score<-rowMeans(a[,paste0(genes,"_logCPM")])
stopifnot(!anyDuplicated(a[,c("sample_id","final_analysis")]))
pre<-a[,c("sample_id","final_analysis","n_cells","state_score")]
names(pre)[3:4]<-c("pre_n","pre_score")
post<-pre; names(post)[3:4]<-c("post_n","post_score")
st<-merge(pairs,pre,by.x="pre_sample_id",by.y="sample_id",all=FALSE)
st<-merge(st,post,by.x=c("post_sample_id","final_analysis"),by.y=c("sample_id","final_analysis"),all=FALSE)
st<-st[st$pre_n>=20 & st$post_n>=20,]
st$delta_state_score<-st$post_score-st$pre_score
sp<-aggregate(delta_state_score~patient+remission+final_analysis,data=st,FUN=mean)
exactp<-function(y,g){
 ok<-is.finite(y)&g%in%c(0,1); y<-y[ok]; g<-g[ok]; n<-length(y); k<-sum(g==1)
 if(n<4||k==0||k==n)return(NA_real_)
 obs<-mean(y[g==1])-mean(y[g==0]); cmb<-combn(n,k)
 dif<-apply(cmb,2,function(ix)mean(y[ix])-mean(y[-ix]))
 mean(abs(dif)>=abs(obs)-1e-12)
}
fitone<-function(state){
 d<-sp[sp$final_analysis==state,]; d$g<-as.integer(d$remission=="Remission")
 if(nrow(d)<4||length(unique(d$g))<2)return(data.frame(final_analysis=state,n=nrow(d),
   remission=sum(d$g),nonremission=sum(1-d$g),mean_remission=NA,mean_nonremission=NA,
   estimate=NA,lower=NA,upper=NA,p_OLS=NA,p_HC3=NA,p_exact_permutation=NA))
 fit<-lm(d$delta_state_score~d$g); X<-model.matrix(fit); e<-resid(fit); h<-hatvalues(fit)
 br<-solve(crossprod(X)); V<-br%*%crossprod(X,X*(e^2/(1-h)^2))%*%br
 b<-coef(fit)[2]; se<-summary(fit)$coef[2,2]; seh<-sqrt(V[2,2]); df<-fit$df.residual; cr<-qt(.975,df)
 data.frame(final_analysis=state,n=nrow(d),remission=sum(d$g),nonremission=sum(1-d$g),
  mean_remission=mean(d$delta_state_score[d$g==1]),mean_nonremission=mean(d$delta_state_score[d$g==0]),
  estimate=b,lower=b-cr*se,upper=b+cr*se,p_OLS=summary(fit)$coef[2,4],
  p_HC3=2*pt(abs(b/seh),df,lower.tail=FALSE),p_exact_permutation=exactp(d$delta_state_score,d$g))
}
states<-sort(unique(a$final_analysis))
res<-do.call(rbind,lapply(states,fitone))
res$BH_fixed_state_family_n15<-p.adjust(res$p_exact_permutation,"BH",n=length(states))
raw<-read.delim(raw_path,check.names=FALSE,stringsAsFactors=FALSE)
raw<-raw[,c("branch","estimate","lower","upper","p_exact_permutation","BH_fixed_state_family_n15")]
names(raw)<-c("final_analysis","raw_estimate","raw_lower","raw_upper","raw_p_exact","raw_BH_n15")
res<-merge(res,raw,by="final_analysis",all.x=TRUE)
res$estimate_change_after_decontx<-res$estimate-res$raw_estimate
write.table(st,file.path(outdir,"decontx_state_site_pairs.tsv"),sep="	",quote=FALSE,row.names=FALSE)
write.table(sp,file.path(outdir,"decontx_state_patient_deltas.tsv"),sep="	",quote=FALSE,row.names=FALSE)
write.table(res,file.path(outdir,"decontx_fixed15_state_models.tsv"),sep="	",quote=FALSE,row.names=FALSE)
targets<-res[res$final_analysis%in%c("Non ileal CT colonocyte","Non ileal TA","Non ileal LGR5pos stem"),]
decision<-list(status="COMPLETE",role="post-result ambient-RNA sensitivity; not a replacement validation",
 fixed_family_n=length(states),evaluable_states=sum(is.finite(res$p_exact_permutation)),
 threshold=20,targets=targets,
 interpretation_rule=list(
  all_three_q_lt_0_05="ambient-robust support for CT, TA and LGR5 state-family extension",
  CT_only_q_lt_0_05="retain CT within-state result; remove cross-lineage claim",
  CT_not_q_lt_0_05="headline CT-state result is ambient-sensitive"),
 created_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
jsonlite::write_json(decision,file.path(outdir,"analysis_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
cat("DECONTX_FIXED15_COMPLETE")
