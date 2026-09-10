options(stringsAsFactors=FALSE,width=160)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(limma)})
root<-'/root/projects/UC_Treatment_Recovery';setwd(root);out<-'provenance/reviews'
o<-readRDS('runs/001A_qc/discovery_preprocessed_v1.rds');y<-o$delta;p<-o$pairs
sets<-readRDS('inputs/annotation/gene_sets_001B.rds')$Hallmark
idx<-lapply(sets,function(g)which(rownames(y)%in%g));idx<-idx[lengths(idx)>=15&lengths(idx)<=500]
H<-t(vapply(idx,function(ii)colMeans(y[ii,,drop=FALSE]),numeric(ncol(y))));colnames(H)<-p$subject
hd<-fread('runs/001B2_calibration/Hallmark_patient_delta.tsv');stopifnot(identical(hd$subject,p$subject))
herror<-max(abs(t(H)-as.matrix(hd[,names(idx),with=FALSE])))
arm<-as.numeric(p$arm=='golimumab');response<-as.numeric(p$response=='Yes')
arm_idx<-list(Placebo=which(arm==0),golimumab=which(arm==1))
manual_ci<-function(e,se,df){e<-unname(e);se<-unname(se);df<-unname(df);c(estimate=e,se=se,lower=e-qt(.975,df)*se,upper=e+qt(.975,df)*se,p=2*pt(-abs(e/se),df))}
ols_ci<-function(z,M,j){b<-qr.solve(M,z);r<-z-M%*%b;df<-length(z)-ncol(M);v<-sum(r^2)/df*solve(crossprod(M));manual_ci(b[j],sqrt(v[j,j]),df)}
effects<-fread('runs/001B2_calibration/Hallmark_patient_effects_CI.tsv');audit<-list();rows<-list()
for(s in rownames(H)){
 z<-H[s,];calc<-list()
 for(a in names(arm_idx)){
  ii<-arm_idx[[a]];za<-z[ii]
  calc[[paste0(a,'_within_change')]]<-manual_ci(mean(za),sd(za)/sqrt(length(za)),length(za)-1)
  calc[[paste0(a,'_response_difference')]]<-ols_ci(za,cbind(1,response[ii]),2)
 }
 x<-z[arm==1];v<-z[arm==0];se2<-var(x)/length(x)+var(v)/length(v);df<-se2^2/((var(x)/length(x))^2/(length(x)-1)+(var(v)/length(v))^2/(length(v)-1))
 calc[['GLM_minus_Placebo']]<-manual_ci(mean(x)-mean(v),sqrt(se2),df)
 calc[['pooled_response_arm_adjusted']]<-ols_ci(z,cbind(1,arm,response),3)
 calc[['arm_response_interaction']]<-ols_ci(z,cbind(1,arm,response,arm*response),4)
 for(co in names(calc)){
  tt<-effects[effects$set==s&effects$contrast==co];stopifnot(nrow(tt)==1)
  observed<-as.numeric(unlist(tt[,names(calc[[co]]),with=FALSE]));d<-max(abs(observed-calc[[co]]))
  audit[[paste(s,co)]]<-d
  rows[[length(rows)+1]]<-data.table(set=s,contrast=co,t(as.matrix(calc[[co]])))
 }
}
re<-rbindlist(rows);fwrite(re,file.path(out,'calibration_effects_independent_v1.tsv'),sep='\t')
manual_bh<-function(p){ii<-order(p);v<-pmin(1,rev(cummin(rev(p[ii]*length(p)/seq_along(p)))));q<-numeric(length(p));q[ii]<-v;q}
effects_bh_error<-max(vapply(unique(effects$contrast),function(co){t<-effects[effects$contrast==co];max(abs(manual_bh(t$p)-t$FDR))},numeric(1)))
# Independently reconstruct pooled model from source fields and evaluate all genes.
bm<-fread('provenance/reviews/metadata_independent_v1_gsm_barcode.tsv')
bl<-bm$inferred_cel_barcode[match(p$baseline_gsm,bm$gsm)];po<-bm$inferred_cel_barcode[match(p$post_gsm,bm$gsm)];lv<-sort(unique(c(bl,po)))
B<-sapply(lv[-1],function(k)as.integer(po==k)-as.integer(bl==k))
pooled<-list();method_repeats<-list()
all<-fread('runs/001B2_calibration/all_Hallmark_tests.tsv');ms<-fread('runs/001B2_calibration/method_sensitivity_counts.tsv');method_count_mismatches<-list();path_bh_errors<-c();five_bh_errors<-c();mixed_bh_errors<-c()
for(br in c('primary_paired','barcode_sensitivity')){
 bdir<-file.path('runs/001B2_calibration',br);W<-cbind(1,arm,response);if(br=='barcode_sensitivity')W<-cbind(W,B)
 d<-fread(file.path(bdir,'pooled_design.tsv'));stopifnot(identical(d$subject,p$subject));design_error<-max(abs(W-as.matrix(d[,-1])))
 beta<-qr.solve(W,t(y));res<-t(y)-W%*%beta;df<-nrow(W)-ncol(W);sigma<-sqrt(colSums(res^2)/df);unscaled<-sqrt(solve(crossprod(W))[3,3])
 fit<-readRDS(file.path(bdir,'pooled_response_fit.rds'));tt<-fread(file.path(bdir,'pooled_response_genes.tsv'));j<-match(rownames(y),tt$gene)
 tval<-beta[3,]/(unscaled*sqrt(fit$s2.post));pv<-2*pt(-abs(tval),df=fit$df.total);half<-qt(.975,fit$df.total)*unscaled*sqrt(fit$s2.post)
 pooled[[br]]<-list(design_max_error=design_error,rank=qr(W)$rank,coefficient_max_error=max(abs(beta[3,]-tt$response_contrast_log2[j])),sigma_max_error=max(abs(sigma-fit$sigma)),P_max_error=max(abs(pv-tt$P.Value[j])),CI_max_error=max(abs(c(beta[3,]-half-tt$CI.L[j],beta[3,]+half-tt$CI.R[j]))),BH_max_error=max(abs(manual_bh(tt$P.Value)-tt$adj.P.Val)),fdr05=sum(tt$adj.P.Val<.05))
 for(met in c('camera_fixed001','camera_estimated','camera_ranks001','fry')){
  block<-all[all$branch==br&all$method==met];five_bh_errors<-c(five_bh_errors,max(abs(manual_bh(block$PValue)-block$FDR_across_five_contrasts)))
  for(co in unique(block$contrast)){
   t<-block[block$contrast==co];path_bh_errors<-c(path_bh_errors,max(abs(manual_bh(t$PValue)-t$FDR)))
   if(met=='fry')mixed_bh_errors<-c(mixed_bh_errors,max(abs(manual_bh(t$PValue.Mixed)-t$FDR.Mixed)))
   expected<-ms[ms$branch==br&ms$method==met&ms$contrast==co]
   if(expected$sets!=nrow(t)||expected$FDR05!=sum(t$FDR<.05))method_count_mismatches[[paste(br,met,co)]]<-TRUE
  }
  # Repeat all 50 Hallmark sets for two pre-specified contrasts per branch/method.
  for(co in c('GLM_change','pooled_response_arm_adjusted')){
   if(co=='GLM_change'){
    des<-as.matrix(fread(file.path('runs/001B_discovery',br,'main_design.tsv'))[,-1]);cv<-setNames(rep(0,ncol(des)),colnames(des));cv['GLM']<-1
   }else{des<-W;cv<-rep(0,ncol(des));cv[3]<-1}
   if(met=='fry')z<-fry(y,idx,design=des,contrast=cv,standardize='residual.sd',sort='none') else z<-camera(y,idx,design=des,contrast=cv,inter.gene.cor=if(met=='camera_estimated')NA else .01,use.ranks=met=='camera_ranks001',sort=FALSE)
   t<-all[all$branch==br&all$method==met&all$contrast==co];j<-match(rownames(z),t$set)
   method_repeats[[paste(br,met,co)]]<-list(P_max_error=max(abs(z$PValue-t$PValue[j])),FDR_max_error=max(abs(z$FDR-t$FDR[j])),direction_identical=identical(as.character(z$Direction),t$Direction[j]),FDR05=sum(z$FDR<.05),correlation_max_error=if(met=='camera_estimated')max(abs(z$Correlation-t$Correlation[j]))else NULL,mixed_P_max_error=if(met=='fry')max(abs(z$PValue.Mixed-t$PValue.Mixed[j]))else NULL)
  }
 }
}
# Recompute raw-probe correlations from centred dot products, independently of cor().
x<-o$probe_expression;m<-o$metadata;base<-which(m$time=='Week 0'&m$disease=='Ulcerative Colitis (UC)');post<-which(m$time=='Week 6'&m$disease=='Ulcerative Colitis (UC)')
xb<-sweep(x[,base],2,colMeans(x[,base]),'-');xp<-sweep(x[,post],2,colMeans(x[,post]),'-')
cc<-crossprod(xb,xp)/outer(sqrt(colSums(xb^2)),sqrt(colSums(xp^2)));rownames(cc)<-m$gsm[base];colnames(cc)<-m$gsm[post]
ic<-fread('runs/001C_programs/identity_qc/all_pre_post_probe_correlations.tsv');ij<-cbind(match(ic$baseline_gsm,rownames(cc)),match(ic$post_gsm,colnames(cc)));corerror<-max(abs(cc[ij]-ic$probe_pearson))
known<-ic[ic$source_ID_match==TRUE];marker<-fread('runs/001C_programs/identity_qc/sex_marker_expression.tsv');marker_error<-max(vapply(c('XIST','RPS4Y1','KDM5D'),function(g)max(abs(o$gene_expression[g,marker$gsm]-marker[[g]])),numeric(1)))
rank_mismatches<-list();rank_all<-c();rank_arm<-c();rank_center<-c()
for(j in seq_along(post)){
 post_index<-post[j];r<-cc[,j];sa<-m$treatment[base]==m$treatment[post_index];sc<-sub('-.*','',m$subject[base])==sub('-.*','',m$subject[post_index]);tab<-ic[ic$post_gsm==m$gsm[post_index]];tab<-tab[match(m$gsm[base],tab$baseline_gsm)]
 ranks<-1+vapply(r,function(v)sum(r>v),integer(1));ar<-rep(NA_real_,length(r));cr<-ar
 ar[sa]<-1+vapply(r[sa],function(v)sum(r[sa]>v),integer(1));cr[sa&sc]<-1+vapply(r[sa&sc],function(v)sum(r[sa&sc]>v),integer(1))
 if(!identical(as.numeric(ranks),as.numeric(tab$correlation_rank_all))||!isTRUE(all.equal(ar,tab$correlation_rank_same_arm,check.attributes=FALSE))||!isTRUE(all.equal(cr,tab$correlation_rank_same_center_arm,check.attributes=FALSE)))rank_mismatches[[m$gsm[post_index]]]<-TRUE
 kk<-which(m$subject[base]==m$subject[post_index]);if(length(kk)){rank_all<-c(rank_all,ranks[kk]);rank_arm<-c(rank_arm,ar[kk]);rank_center<-c(rank_center,cr[kk])}
}
is<-fromJSON('runs/001C_programs/identity_qc/identity_qc_summary.json')
identity<-list(all_pairwise_correlations_checked=nrow(ic),correlation_max_error=corerror,marker_max_error=marker_error,rank_mismatches=rank_mismatches,documented_pairs=length(rank_all),top_all=sum(rank_all==1),top_same_arm=sum(rank_arm==1),top_same_center_arm=sum(rank_center==1),declared_new_pairs=is$additional_pairs_accepted,known_vs_other_cor_ranges=list(known=range(ic$probe_pearson[ic$source_ID_match]),other=range(ic$probe_pearson[!ic$source_ID_match])))
report<-list(status='COMPLETE',current_directory_only='runs/001B2_calibration; deprecated export_v1 never read',Hallmark_patient_delta_max_error=herror,effects_checked=length(audit),effect_CI_max_error=max(unlist(audit)),patient_effect_BH_max_error=effects_bh_error,pooled_gene_model_checks=pooled,method_repeat_checks=method_repeats,method_summary_count_mismatches=method_count_mismatches,all_test_BH_max_error=max(path_bh_errors),five_contrast_BH_max_error=max(five_bh_errors),fry_mixed_BH_max_error=max(mixed_bh_errors),identity=identity,session=capture.output(sessionInfo()))
write_json(report,file.path(out,'calibration_identity_numeric_v1.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat(toJSON(list(score_error=herror,effects_checked=length(audit),effects_max_error=max(unlist(audit)),pooled=pooled,method_repeat_max_P_error=max(vapply(method_repeats,`[[`,numeric(1),'P_max_error')),identity=identity),pretty=TRUE,auto_unbox=TRUE,digits=16),'\n')
