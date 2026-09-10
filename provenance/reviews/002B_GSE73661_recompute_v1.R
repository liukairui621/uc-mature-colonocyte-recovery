options(stringsAsFactors=FALSE,width=180)
suppressPackageStartupMessages(library(jsonlite))
setwd("/root/projects/UC_Treatment_Recovery");out<-"provenance/reviews";prod<-"runs/002B_GSE73661_support"
plan<-read_json("planning/support_plan_002B.frozen.json",simplifyVector=TRUE);lock<-read_json("signature_lock.json",simplifyVector=TRUE)
stopifnot(plan$status=="FROZEN",plan$role=="descriptive_supportive_only",plan$n_confirmatory_primary_tests==0)
con<-gzfile("inputs/support/GSE73661_series_matrix.txt.gz","rt");lines<-readLines(con,warn=FALSE);close(con)
beg<-which(lines=="!series_matrix_table_begin");end<-which(lines=="!series_matrix_table_end");stopifnot(length(beg)==1,length(end)==1,end>beg)
tab<-read.table(text=paste(lines[(beg+1):(end-1)],collapse="\n"),header=TRUE,sep="\t",quote='"',comment.char="",check.names=FALSE,colClasses=c("character",rep("numeric",178)))
E<-as.matrix(tab[,-1]);rownames(E)<-tab[,1];rm(tab)
stopifnot(identical(dim(E),c(33252L,178L)),all(is.finite(E)),!anyDuplicated(rownames(E)),!anyDuplicated(colnames(E)))
q<-quantile(E,c(0,.01,.25,.5,.75,.99,1));stopifnot(q[6]<=40,q[1]>=-20)
header<-lapply(lines[seq_len(beg-1)],function(z)scan(text=z,what="",sep="\t",quote='"',quiet=TRUE));rm(lines)
charrows<-header[vapply(header,function(z)length(z)>0&&z[1]=="!Sample_characteristics_ch1",logical(1))]
field<-function(prefix){z<-charrows[vapply(charrows,function(z)startsWith(z[2],paste0(prefix,": ")),logical(1))];stopifnot(length(z)==1);substring(z[[1]][-1],nchar(prefix)+3)}
gsm<-header[vapply(header,function(z)length(z)>0&&z[1]=="!Sample_geo_accession",logical(1))][[1]][-1]
md<-data.frame(gsm=gsm,subject=field("study individual number"),visit=field("week (w)"),treatment=field("induction therapy_maintenance therapy"),mayo=field("mayo endoscopic subscore"))
stopifnot(identical(md$gsm,colnames(E)));ifx<-md[md$treatment=="IFX",];subjects<-unique(ifx$subject);subjects<-subjects[order(as.integer(subjects))]
pairs<-do.call(rbind,lapply(subjects,function(s){z<-ifx[ifx$subject==s,];stopifnot(sum(z$visit=="W0")==1,sum(z$visit=="W4_W6")==1);bm<-as.integer(z$mayo[z$visit=="W0"]);pm<-as.integer(z$mayo[z$visit=="W4_W6"]);data.frame(subject=s,baseline_gsm=z$gsm[z$visit=="W0"],followup_gsm=z$gsm[z$visit=="W4_W6"],healing=ifelse(pm<=1,"Yes","No"),baseline_endoscopic=bm,post_endoscopic=pm)}))
stopifnot(nrow(pairs)==23,sum(pairs$healing=="Yes")==8)
frozen<-read.delim(plan$pair_file,check.names=FALSE)
stopifnot(identical(pairs$subject,as.character(frozen$subject)),identical(pairs$baseline_gsm,frozen$baseline_gsm),identical(pairs$followup_gsm,frozen$followup_gsm),identical(pairs$healing,frozen$endoscopic_healing),identical(pairs$baseline_endoscopic,frozen$baseline_mayo_endoscopic),identical(pairs$post_endoscopic,frozen$followup_mayo_endoscopic))
map<-read.delim("runs/001D_annotation/GPL6244_canonical_probe_mapping.tsv",check.names=FALSE,na.strings="");stopifnot(!anyDuplicated(map$probe))
i<-match(rownames(E),map$probe);keep<-!is.na(i)&map$retained[i]==1&!is.na(map$canonical_symbol[i]);keep[is.na(keep)]<-FALSE;sym<-map$canonical_symbol[i[keep]]
sumG<-rowsum(E[keep,,drop=FALSE],group=sym,reorder=TRUE);num<-table(sym);G<-sumG/as.numeric(num[rownames(sumG)]);G<-G[order(rownames(G)),,drop=FALSE]
cand<-lock$candidate$id;infl<-"HALLMARK_INFLAMMATORY_RESPONSE";refs<-lock$secondary_reference_programs;sets<-list();sets[[cand]]<-plan$candidate_genes;sets[[infl]]<-plan$inflammation_genes
stopifnot(identical(sets[[cand]],lock$candidate$genes),identical(sets[[infl]],lock$inflammation_covariate$genes),length(sets[[cand]])==6,length(sets[[infl]])==194,all(unlist(sets)%in%rownames(G)))
mem<-read.delim("runs/001D_annotation/program_membership_canonical.tsv",check.names=FALSE);covs<-list()
for(s in c(cand,infl,refs)){src<-unique(mem$gene[mem$program==s]);hit<-intersect(src,rownames(G));covs[[s]]<-data.frame(program=s,source_n=length(src),measured_n=length(hit),coverage=length(hit)/length(src));if(s%in%refs)sets[[s]]<-hit}
W<-matrix(0,nrow=length(sets),ncol=nrow(G),dimnames=list(names(sets),rownames(G)));for(s in names(sets))W[s,sets[[s]]]<-1/length(sets[[s]])
S<-W%*%G;BL<-S[,pairs$baseline_gsm,drop=FALSE];PO<-S[,pairs$followup_gsm,drop=FALSE];D<-PO-BL;colnames(BL)<-colnames(PO)<-colnames(D)<-pairs$subject
stored<-readRDS(file.path(prod,"support_scores.rds"))
precomp<-list(deposited_log2_probe=max(abs(E-stored$probe_expression_log2)),all_gene_means=max(abs(G-stored$gene_expression[rownames(G),colnames(G)])),all_program_scores=max(abs(S-stored$sample_scores[rownames(S),colnames(S)])),all_baseline=max(abs(BL-stored$baseline[rownames(BL),colnames(BL)])),all_deltas=max(abs(D-stored$delta[rownames(D),colnames(D)])))
stopifnot(max(unlist(precomp))<1e-10)
makeX<-function(s,model){X<-cbind(Intercept=1,endoscopic_healingYes=as.numeric(pairs$healing=="Yes"));if(model!="parallel_no_baseline")X<-cbind(X,baseline=as.numeric(BL[s,]));if(model!="parallel_no_inflammation")X<-cbind(X,inflammation=as.numeric(D[infl,]));X}
svdmodel<-function(X,y){
 z<-svd(X);rank<-sum(z$d>max(z$d)*max(dim(X))*.Machine$double.eps);stopifnot(rank==ncol(X));beta<-as.vector(z$v%*%((crossprod(z$u,y))/z$d));names(beta)<-colnames(X)
 bread<-tcrossprod(sweep(z$v,2,z$d,"/"));dimnames(bread)<-list(colnames(X),colnames(X));res<-as.vector(y-X%*%beta);df<-length(y)-rank;sigma<-sqrt(sum(res^2)/df);se<-sqrt(diag(bread))*sigma;h<-rowSums(z$u^2);A<-bread%*%t(X);robse<-sqrt(rowSums(sweep(A,2,res/(1-h),"*")^2));crit<-qt(.975,df);term<-"endoscopic_healingYes"
 vifs<-diag(solve(cor(X[,-1,drop=FALSE])));k<-kappa(cbind(1,scale(X[,-1,drop=FALSE])),exact=TRUE)
 row<-function(sem,method)data.frame(method=method,estimate=unname(beta[term]),se=unname(sem[term]),lower=unname(beta[term]-crit*sem[term]),upper=unname(beta[term]+crit*sem[term]),p=unname(2*pt(-abs(beta[term]/sem[term]),df)),n=length(y),df=df,rank=rank,standardized_condition=k,max_VIF=max(vifs),max_leverage=max(h))
 list(OLS=row(se,"OLS"),HC3=row(robse,"HC3"),beta=beta,res=res,h=h,sigma=sigma,X=X,VIF=vifs)
}
fits<-list();effects<-list()
for(label in c("adjusted_support","parallel_no_baseline","parallel_no_inflammation")){f<-svdmodel(makeX(cand,label),as.numeric(D[cand,]));fits[[label]]<-f;effects[[label]]<-cbind(program=cand,model=label,rbind(f$OLS,f$HC3))}
effects<-do.call(rbind,effects)
compare<-function(own,file,keys,ncols){obs<-read.delim(file.path(prod,file),check.names=FALSE);key<-function(z)apply(z[,keys,drop=FALSE],1,paste,collapse="|");ka<-key(own);kb<-key(obs);stopifnot(!anyDuplicated(ka),!anyDuplicated(kb),setequal(ka,kb));mx<-max(abs(as.matrix(own[,ncols,drop=FALSE])-as.matrix(obs[match(ka,kb),ncols,drop=FALSE])));stopifnot(mx<1e-8);list(file=file,n_rows=nrow(own),max_numeric_difference=mx)}
numcols<-setdiff(names(fits[[1]]$OLS),"method");checks<-list(compare(effects,"candidate_support_models.tsv",c("program","model","method"),numcols))
sec<-do.call(rbind,lapply(refs,function(s){stopifnot(covs[[s]]$coverage>=.8);f<-svdmodel(makeX(s,"adjusted_support"),as.numeric(D[s,]));cbind(program=s,model="reference_descriptive",f$OLS)}))
o<-order(sec$p);vals<-pmin(1,rev(cummin(rev(sec$p[o]*3/seq_along(o)))));sec$BH_across_three_descriptive<-NA_real_;sec$BH_across_three_descriptive[o]<-vals
checks[[length(checks)+1]]<-compare(sec,"reference_descriptive_models.tsv",c("program","model","method"),c(numcols,"BH_across_three_descriptive"))
f<-fits$adjusted_support;X<-f$X;y<-as.numeric(D[cand,]);lo<-do.call(rbind,lapply(seq_along(y),function(i){z<-svdmodel(X[-i,,drop=FALSE],y[-i]);data.frame(omitted_subject=pairs$subject[i],estimate=z$OLS$estimate,p=z$OLS$p)}));checks[[length(checks)+1]]<-compare(lo,"leave_one_patient_out.tsv","omitted_subject",c("estimate","p"))
df<-f$OLS$df;mse<-f$sigma^2;extsigma<-sqrt((sum(f$res^2)-f$res^2/(1-f$h))/(df-1));diag<-data.frame(subject=pairs$subject,leverage=f$h,cooks_distance=(f$res^2/(ncol(X)*mse))*f$h/(1-f$h)^2,residual=f$res,studentized_residual=f$res/(extsigma*sqrt(1-f$h)))
checks[[length(checks)+1]]<-compare(diag,"patient_diagnostics.tsv","subject",c("leverage","cooks_distance","residual","studentized_residual"))
plotdata<-data.frame(subject=pairs$subject,baseline_candidate=as.numeric(BL[cand,]),post_candidate=as.numeric(PO[cand,]),delta_candidate=as.numeric(D[cand,]),delta_inflammation=as.numeric(D[infl,]),baseline_endoscopic=pairs$baseline_endoscopic,post_endoscopic=pairs$post_endoscopic)
checks[[length(checks)+1]]<-compare(plotdata,"paired_score_plot_source.tsv","subject",setdiff(names(plotdata),"subject"))
response<-X[,"endoscopic_healingYes"];Z<-X[,colnames(X)!="endoscopic_healingYes",drop=FALSE];z<-svd(Z);b<-z$v%*%((crossprod(z$u,response))/z$d);rssr<-sum((response-Z%*%b)^2);vifr<-sum((response-mean(response))^2)/rssr
stopifnot(abs(f$OLS$se-f$sigma/sqrt(rssr))<1e-10)
summary<-read_json(file.path(prod,"support_summary.json"),simplifyVector=TRUE);stopifnot(abs(f$sigma-summary$residual_sigma)<1e-10,abs(sd(y)-summary$candidate_delta_SD)<1e-10,max(abs(as.numeric(q)-unlist(summary$source_quantiles)))<1e-10)
stopifnot(f$OLS$n==23,f$OLS$rank==4,f$OLS$df==19,summary$GSE23597_primary_support_remains_FALSE,!summary$can_replace_GSE23597,!summary$pooled_P_computed)
write.table(effects,file.path(out,"002B_GSE73661_recomputed_candidate_v1.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(sec,file.path(out,"002B_GSE73661_recomputed_reference_v1.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(diag,file.path(out,"002B_GSE73661_recomputed_diagnostics_v1.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(lo,file.path(out,"002B_GSE73661_recomputed_LOPO_v1.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
result<-list(status="INDEPENDENT_RECOMPUTATION_PASSED",role="descriptive_supportive_only",raw_dimensions=dim(E),source_quantiles=as.list(q),log_transform_applied=FALSE,normalization_applied=FALSE,retained_probe_count=sum(keep),canonical_gene_count=nrow(G),IFX_pairs=nrow(pairs),healing=sum(pairs$healing=="Yes"),nonhealing=sum(pairs$healing=="No"),coverage=do.call(rbind,covs),scored_common_inflammation=length(sets[[infl]]),preprocessing_max_differences=precomp,model_table_checks=checks,candidate_models=effects,reference_models=sec,adjusted_residual_sigma=f$sigma,observed_candidate_delta_SD=sd(y),response_residualized_SS=unname(rssr),response_VIF=unname(vifr),reconstructed_SE=f$sigma/sqrt(rssr),LOPO_estimate_range=range(lo$estimate),LOPO_p_range=range(lo$p),LOPO_all_full_rank=TRUE,maximum_leverage_patient=diag[which.max(diag$leverage),],maximum_Cooks_patient=diag[which.max(diag$cooks_distance),],high_leverage_descriptive=diag[diag$leverage>2*ncol(X)/nrow(X),],GSE23597_primary_support_remains_FALSE=TRUE,confirmatory_primary_test_count=0,methods="Independent raw-table parser and raw-characteristic IFX pairing; deposited RMA values retained; canonical rowsum/probe-count aggregation; explicit membership weight matrix; SVD OLS/covariance, direct sandwich-column HC3, manual BH, analytical patient diagnostics and all LOPO. No limma avereps or lm call in independent recomputation.")
write_json(result,file.path(out,"002B_GSE73661_numeric_independent_v1.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
print(effects);print(sec);print(diag[order(-diag$leverage),][1:5,]);print(result[c("adjusted_residual_sigma","response_VIF","reconstructed_SE","LOPO_estimate_range","LOPO_p_range")]);print(precomp);cat("GSE73661_INDEPENDENT_RECOMPUTATION_PASSED\n")
