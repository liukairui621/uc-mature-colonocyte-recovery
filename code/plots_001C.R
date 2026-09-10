options(stringsAsFactors=FALSE);set.seed(20260910)
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
root<-"/root/projects/UC_Treatment_Recovery";setwd(root);figdir<-"runs/001C_figures";dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
source("code/figure_theme.R")
a<-fread("runs/001B2_calibration/Hallmark_patient_effects_CI.tsv",sep="\t")
labs<-c(HALLMARK_INTERFERON_GAMMA_RESPONSE="IFN-gamma",HALLMARK_INTERFERON_ALPHA_RESPONSE="IFN-alpha",HALLMARK_IL6_JAK_STAT3_SIGNALING="IL6 / JAK / STAT3",HALLMARK_INFLAMMATORY_RESPONSE="Inflammatory response",HALLMARK_OXIDATIVE_PHOSPHORYLATION="Oxidative phosphorylation")
w<-a[set%in%names(labs)&contrast%in%c("golimumab_within_change","Placebo_within_change")]
w[,arm:=ifelse(contrast=="golimumab_within_change","golimumab","Placebo")]
w[,label:=factor(unname(labs[set]),levels=rev(unname(labs)))]
fwrite(w,file.path(figdir,"Revision01_WithinArm_CI_source.tsv"),sep="\t")
po<-position_dodge(width=.45)
p<-ggplot(w,aes(x=estimate,y=label,color=arm))+geom_vline(xintercept=0,color="grey75",linewidth=.4)+geom_errorbar(aes(xmin=lower,xmax=upper),width=.16,position=po,orientation="y",linewidth=.6)+geom_point(position=po,size=2.3)+scale_color_manual(values=project_palette[c("golimumab","Placebo")],labels=c(golimumab="Golimumab (42 patients)",Placebo="Placebo (23 patients)"))+labs(x="Mean paired change in log2 expression (95% CI)",y=NULL,color=NULL,title="Five Hallmark programs highlighted in the review",subtitle="Patient-level means; group-specific significance is not a treatment comparison")+theme_project()
save_panel(p,"Revision01_WithinArm_CI",8,4.4)
d<-a[set%in%names(labs)&contrast=="GLM_minus_Placebo"];d[,label:=factor(unname(labs[set]),levels=rev(unname(labs)))]
fwrite(d,file.path(figdir,"Revision02_ArmDifference_CI_source.tsv"),sep="\t")
p<-ggplot(d,aes(x=estimate,y=label))+geom_vline(xintercept=0,color="grey75",linewidth=.4)+geom_errorbar(aes(xmin=lower,xmax=upper),width=.16,orientation="y",linewidth=.6)+geom_point(size=2.3,color="#2166AC")+labs(x="Golimumab minus placebo: difference in log2 change (95% CI)",y=NULL,title="Direct comparison between treatment arms",subtitle="Welch intervals; intervals spanning zero do not establish equivalence")+theme_project()
save_panel(p,"Revision02_ArmDifference_CI",8,4.3)
ef<-fread("runs/001C2_refinement/refinement_response_effects_CI.tsv",sep="\t")
ml<-c(common_baseline="Arm + baseline program",common_baseline_inflammation="+ inflammation change",common_baseline_inflammation_epithelial="+ epithelial marker proxy",common_baseline_inflammation_barcode="+ inferred barcode",common_baseline_inflammation_MARS_TableS1_11="+ MARS printed-table components",common_baseline_inflammation_MARS_v74_11="+ MARS v7.4 components",candidate_conditioned_on_CT_program="+ CT-colonocyte change and baseline")
z<-ef[program=="MATURE_COLONOCYTE6_EXPLORATORY"];z[,model_label:=factor(unname(ml[model]),levels=rev(unname(ml)))]
fwrite(z,file.path(figdir,"Exploratory03_CandidateConditional_CI_source.tsv"),sep="\t")
p<-ggplot(z,aes(x=estimate,y=model_label))+geom_vline(xintercept=0,color="grey75",linewidth=.4)+geom_errorbar(aes(xmin=lower,xmax=upper),width=.14,orientation="y",linewidth=.6)+geom_point(size=2.3,color="#2166AC")+labs(x="Responder-associated difference in candidate log2 change (95% CI)",y=NULL,title="Six-gene candidate: conditional associations",subtitle="Each extension starts from the arm + baseline + inflammation model;\nMARS comparisons use component means, not the original MARS classifier")+theme_project()
save_panel(p,"Exploratory03_CandidateConditional_CI",9,5.1)
o<-readRDS("runs/001C_programs/program_scores.rds");s<-"MATURE_COLONOCYTE6_EXPLORATORY";ct<-"UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy"
df<-as.data.table(o$pairs);df[,candidate_delta:=o$delta[s,]]
long<-rbind(data.table(df,reference="Inflammatory response",reference_delta=as.numeric(o$delta["HALLMARK_INFLAMMATORY_RESPONSE",])),data.table(df,reference="Published CT-colonocyte markers (18/20)",reference_delta=as.numeric(o$delta[ct,])))
fwrite(long,file.path(figdir,"Exploratory04_CandidateRelations_source.tsv"),sep="\t")
p<-ggplot(long,aes(x=reference_delta,y=candidate_delta,color=response,shape=arm))+geom_hline(yintercept=0,color="grey85",linewidth=.35)+geom_vline(xintercept=0,color="grey85",linewidth=.35)+geom_point(size=2.2,alpha=.8)+facet_wrap(~reference,nrow=1,scales="free_x")+scale_color_manual(values=project_palette[c("No","Yes")])+scale_shape_manual(values=c(Placebo=1,golimumab=16))+labs(x="Reference-program change (mean log2)",y="Six-gene candidate change (mean log2)",color="Clinical response",shape="Treatment arm",title="A known mature-epithelial axis, with patient-level variation",subtitle="One point per patient in each panel (n = 65); this is not evidence of a new cell state")+theme_project()
save_panel(p,"Exploratory04_CandidateRelations",10.5,4.6)
capture.output(sessionInfo(),file=file.path(figdir,"sessionInfo.txt"))
cat("FIGURES_COMPLETE\n")
