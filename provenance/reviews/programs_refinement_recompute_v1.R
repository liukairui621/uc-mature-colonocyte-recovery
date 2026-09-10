options(stringsAsFactors=FALSE,width=160)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setwd('/root/projects/UC_Treatment_Recovery');out<-'provenance/reviews'
g<-readRDS('runs/001A_qc/discovery_preprocessed_v1.rds')$gene_expression
o<-readRDS('runs/001C_programs/program_scores.rds');p<-o$pairs
mem<-fread('inputs/annotation/program_membership_001C.tsv');programs<-sort(unique(mem$program))
# Independent matrix-weight reconstruction from original gene-expression RDS.
W<-matrix(0,length(programs),nrow(g),dimnames=list(programs,rownames(g)))
for(s in programs){m<-mem[mem$program==s&mem$gene%in%rownames(g)];W[s,m$gene]<-m$weight/sum(abs(m$weight))}
S<-W%*%g;S<-S[rownames(o$sample_scores),colnames(o$sample_scores)]
BL<-S[,p$baseline_gsm];PO<-S[,p$post_gsm];D<-PO-BL;colnames(BL)<-colnames(PO)<-colnames(D)<-p$subject
score_errors<-list(all_sample=max(abs(S-o$sample_scores)),baseline=max(abs(BL-o$baseline)),post=max(abs(PO-o$post)),delta=max(abs(D-o$delta)))
cand<-'MATURE_COLONOCYTE6_EXPLORATORY';CT<-'UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy';ox<-'HALLMARK_OXIDATIVE_PHOSPHORYLATION'
arm<-as.numeric(p$arm=='golimumab');resp<-as.numeric(p$response=='Yes');mayo<-as.numeric(p$mayo_baseline)
inf<-D['HALLMARK_INFLAMMATORY_RESPONSE',];epi<-D['EPITHELIAL_ABUNDANCE_PROXY4',]
bm<-fread('provenance/reviews/metadata_independent_v1_gsm_barcode.tsv');bb<-bm$inferred_cel_barcode[match(p$baseline_gsm,bm$gsm)];bp<-bm$inferred_cel_barcode[match(p$post_gsm,bm$gsm)];lv<-sort(unique(c(bb,bp)));B<-sapply(lv[-1],function(k)as.integer(bp==k)-as.integer(bb==k))
manual_bh<-function(p){i<-order(p);v<-pmin(1,rev(cummin(rev(p[i]*length(p)/seq_along(p)))));q<-numeric(length(p));q[i]<-v;q}
ols<-function(z,X,j=3){b<-qr.solve(X,z);r<-as.numeric(z-X%*%b);df<-nrow(X)-ncol(X);V<-sum(r^2)/df*solve(crossprod(X));se<-sqrt(V[j,j]);e<-unname(b[j]);half<-qt(.975,df)*se;list(values=c(estimate=e,se=se,lower=e-half,upper=e+half,p=2*pt(-abs(e/se),df)),b=b,residuals=r,df=df,X=X)}
record_check<-function(calc,row,fields=c('estimate','se','lower','upper','p'))max(abs(calc[fields]-as.numeric(unlist(row[,fields,with=FALSE]))))
ctab<-fread('runs/001C_programs/program_response_and_arm_effects_CI.tsv');rtab<-fread('runs/001C2_refinement/refinement_response_effects_CI.tsv');cerrors<-c();rerrors<-c();candidate_rows<-list()
for(i in seq_len(nrow(ctab))){
 s<-ctab$program[i];model<-ctab$model[i];z<-D[s,];base<-BL[s,];j<-3
 X<-switch(model,response_arm_adjusted=cbind(1,arm,resp),response_baseline_adjusted=cbind(1,arm,resp,base,mayo),response_barcode_adjusted=cbind(1,arm,resp,B),arm_response_interaction=cbind(1,arm,resp,arm*resp),arm_baseline_adjusted=cbind(1,arm,base,mayo),response_plus_inflammation=cbind(1,arm,resp,base,mayo,inf),response_plus_epithelial_proxy=cbind(1,arm,resp,base,mayo,epi),response_plus_both=cbind(1,arm,resp,base,mayo,inf,epi),NULL)
 if(model=='arm_response_interaction')j<-4
 if(model=='arm_baseline_adjusted')j<-2
 if(model%in%c('Placebo_response_difference','golimumab_response_difference')){ii<-which(arm==as.numeric(startsWith(model,'golimumab')));z<-z[ii];X<-cbind(1,resp[ii]);j<-2}
 a<-ols(z,X,j);cerrors<-c(cerrors,record_check(a$values,ctab[i]))
 if(s==cand)candidate_rows[[paste('001C',model)]]<-as.list(a$values)
}
for(i in seq_len(nrow(rtab))){
 s<-rtab$program[i];model<-rtab$model[i];z<-D[s,];base<-BL[s,];X<-cbind(1,arm,resp,base)
 if(model!='common_baseline')X<-cbind(X,inf)
 if(model=='common_baseline_inflammation_epithelial')X<-cbind(X,epi)
 if(model=='common_baseline_inflammation_barcode')X<-cbind(X,B)
 if(grepl('MARS_TableS1_11$',model))X<-cbind(X,t(D[grep('^MARS_TableS1_column_',rownames(D)),,drop=FALSE]))
 if(grepl('MARS_v74_11$',model))X<-cbind(X,t(D[grep('^MARS_MSigDB74_',rownames(D)),,drop=FALSE]))
 if(model=='candidate_conditioned_on_CT_program')X<-cbind(X,D[CT,],BL[CT,])
 a<-ols(z,X);rerrors<-c(rerrors,record_check(a$values,rtab[i]))
 if(s==cand)candidate_rows[[paste('001C2',model)]]<-as.list(a$values)
}
# Independent Pearson coefficient and Fisher interval from paired scores.
x<-as.numeric(D[cand,]);z<-as.numeric(D[CT,]);xc<-x-mean(x);zc<-z-mean(z);r<-sum(xc*zc)/sqrt(sum(xc^2)*sum(zc^2));ci<-tanh(atanh(r)+c(-1,1)*qnorm(.975)/sqrt(length(x)-3));pval<-2*pt(-abs(r*sqrt((length(x)-2)/(1-r^2))),df=length(x)-2)
cor_table<-fread('runs/001C_programs/program_overlap_and_delta_correlations.tsv');row<-cor_table[cor_table$program==cand&cor_table$reference==CT]
corcheck<-list(r=r,CI=ci,P=pval,max_error=max(abs(c(r,ci,pval)-c(row$pearson,row$lower,row$upper,row$p))),shared_genes=intersect(mem$gene[mem$program==cand],mem$gene[mem$program==CT]))
# Independent Newton/IRLS binomial MLE with explicitly evaluated log likelihood.
loglik<-function(beta,X,y){eta<-as.numeric(X%*%beta);sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))}
logit_fit<-function(X,y){b<-numeric(ncol(X));b[1]<-qlogis(mean(y));converged<-FALSE
 for(it in 1:100){eta<-as.numeric(X%*%b);mu<-plogis(eta);w<-mu*(1-mu);step<-solve(crossprod(X,X*w),crossprod(X,y-mu));bn<-b+as.numeric(step);if(max(abs(bn-b))<1e-11){b<-bn;converged<-TRUE;break};b<-bn}
 mu<-plogis(as.numeric(X%*%b));V<-solve(crossprod(X,X*(mu*(1-mu))));list(beta=b,se=sqrt(diag(V)),loglik=loglik(b,X,y),converged=converged)}
ltab<-fread('runs/001C_programs/incremental_response_association.tsv');ltab2<-fread('runs/001C2_refinement/common_covariate_incremental_association.tsv');logs<-list()
for(stage in c('001C','001C2')){
 t<-if(stage=='001C')ltab else ltab2
 for(i in seq_len(nrow(t))){s<-t$program[i];spec<-if(stage=='001C')t$model[i]else 'common_covariate';X0<-cbind(1,arm,inf)
  if(stage=='001C2')X0<-cbind(X0,BL[s,]) else if(spec=='arm_inflammation_baseline')X0<-cbind(X0,BL[s,],mayo)
  f0<-logit_fit(X0,resp);f1<-logit_fit(cbind(X0,D[s,]),resp);j<-length(f1$beta);be<-f1$beta[j];se<-f1$se[j];lrt<-pchisq(2*(f1$loglik-f0$loglik),df=1,lower.tail=FALSE);or<-exp(be*sd(D[s,]));lo<-exp((be-qnorm(.975)*se)*sd(D[s,]));hi<-exp((be+qnorm(.975)*se)*sd(D[s,]))
  declared<-if(stage=='001C')c(t$beta_per_log2[i],t$SE[i],t$LRT_p[i],t$OR_per_SD[i],t$OR_SD_lower[i],t$OR_SD_upper[i]) else c(t$beta[i],t$se[i],t$LRT_p[i],t$OR_per_SD[i],t$lower_OR[i],t$upper_OR[i])
  logs[[paste(stage,s,spec)]]<-list(beta=be,se=se,LRT_p=lrt,OR_SD=or,OR_CI=c(lo,hi),all_values_max_error=max(abs(c(be,se,lrt,or,lo,hi)-declared)),converged=f0$converged&&f1$converged)
 }
}
# Check all relevant local BH families without interpreting them as global control.
bhchecks<-list();families<-list()
for(model_key in unique(ctab$model))for(fam in unique(ctab$family)){
 t<-ctab[ctab$model==model_key&ctab$family==fam];if(nrow(t)) {bhchecks[[paste('001C',model_key,fam)]]<-max(abs(manual_bh(t$p)-t$FDR));families[[paste('001C',model_key,fam)]]<-nrow(t)}
}
for(model_key in unique(rtab$model)){t<-rtab[rtab$model==model_key];bhchecks[[paste('001C2',model_key)]]<-max(abs(manual_bh(t$p)-t$FDR_across_declared_programs));families[[paste('001C2',model_key)]]<-nrow(t)}
for(model_key in unique(ltab$model)){t<-ltab[ltab$model==model_key];bhchecks[[paste('001C_logistic',model_key)]]<-max(abs(manual_bh(t$LRT_p)-t$FDR))}
bhchecks[['001C2_logistic']]<-max(abs(manual_bh(ltab2$LRT_p)-ltab2$FDR_across_four_programs))
# HC3 computed from independent QR residuals and direct sandwich products.
htab<-fread('runs/001C2_refinement/HC3_common_model_sensitivity.tsv');hc<-list()
for(i in seq_len(nrow(htab))){s<-htab$program[i];model<-htab$model[i];X<-cbind(1,arm,resp,BL[s,]);if(model=='common_baseline_inflammation')X<-cbind(X,inf)
 a<-ols(D[s,],X);bread<-solve(crossprod(X));h<-rowSums((X%*%bread)*X);var3<-bread%*%crossprod(X, X*as.numeric(a$residuals^2/(1-h)^2))%*%bread;se<-sqrt(var3[3,3]);est<-a$values['estimate'];lo<-unname(est-qt(.975,a$df)*se);hi<-unname(est+qt(.975,a$df)*se);pv<-2*pt(-abs(est/se),df=a$df)
 xs<-cbind(1,scale(X[,-1,drop=FALSE]));sv<-svd(xs,nu=0,nv=0)$d;condition<-max(sv)/min(sv)
 hc[[paste(s,model)]]<-list(estimate=unname(est),HC3_se=se,CI=c(lo,hi),P=unname(pv),condition_standardized=condition,max_error=max(abs(c(est,se,lo,hi,pv,condition,max(h))-c(htab$estimate[i],htab$HC3_se[i],htab$HC3_lower[i],htab$HC3_upper[i],htab$HC3_p[i],htab$scaled_design_condition_number[i],htab$max_leverage[i]))))
}
for(model_key in unique(htab$model)){t<-htab[htab$model==model_key];bhchecks[[paste('HC3',model_key)]]<-max(abs(manual_bh(t$HC3_p)-t$HC3_FDR))}
# Recompute patient deletions for the fixed candidate, preserving all observations.
l1<-fread('runs/001C_programs/leave_one_patient_out.tsv');l2<-fread('runs/001C2_refinement/refinement_leave_one_patient_out.tsv');loo<-list()
for(stage in c('001C','001C2')){tab<-if(stage=='001C')l1 else l2;tab<-tab[tab$program==cand];X<-cbind(1,arm,resp);if(stage=='001C2')X<-cbind(X,BL[cand,],inf);errs<-c();values<-c()
 for(i in seq_len(nrow(p))){be<-ols(D[cand,-i],X[-i,,drop=FALSE])$values['estimate'];values<-c(values,be);decl<-tab$estimate[match(p$subject[i],tab$omitted_subject)];errs<-c(errs,abs(be-decl))}
 loo[[stage]]<-list(max_error=max(errs),estimate_range=range(values),all_same_positive_sign=all(values>0))
}
# Coverage reconstructed from exact membership, no imputation or relabelling.
cov<-fread('runs/001C_programs/program_coverage.tsv');coverage<-lapply(programs,function(s){m<-mem[mem$program==s];data.table(program=s,source=nrow(m),measured=sum(m$gene%in%rownames(g)),fraction=mean(m$gene%in%rownames(g))) });coverage<-rbindlist(coverage)
coverage_order<-match(coverage$program,cov$program);cv<-cov[coverage_order];coverage_match<-all(coverage$source==cv$source_genes)&all(coverage$measured==cv$measured_genes)&max(abs(coverage$fraction-cv$coverage))<1e-12
report<-list(status='COMPLETE',scores_max_errors=score_errors,membership_all_weights_one=all(mem$weight==1),candidate_members=mem$gene[mem$program==cand],n_subjects=nrow(p),n_programs=nrow(D),n_001C_OLS_models=length(cerrors),max_001C_OLS_error=max(cerrors),n_001C2_OLS_models=length(rerrors),max_001C2_OLS_error=max(rerrors),candidate_OLS_recomputed=candidate_rows,candidate_CT_correlation=corcheck,logistic_recomputed=logs,BH_family_sizes=families,BH_max_errors=bhchecks,HC3_recomputed=hc,candidate_leave_one_patient_out=loo,coverage_identical=coverage_match,coverage=coverage,session=capture.output(sessionInfo()))
write_json(report,file.path(out,'programs_refinement_numeric_v1.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat(toJSON(list(scores=score_errors,n_models=c(length(cerrors),length(rerrors)),OLS_max_error=max(c(cerrors,rerrors)),candidate_common=candidate_rows[['001C2 common_baseline_inflammation']],candidate_CT=corcheck,max_logistic_error=max(vapply(logs,`[[`,numeric(1),'all_values_max_error')),BH_max_error=max(unlist(bhchecks)),HC3_max_error=max(vapply(hc,`[[`,numeric(1),'max_error')),coverage_match=coverage_match),pretty=TRUE,auto_unbox=TRUE,digits=16),'\n')
