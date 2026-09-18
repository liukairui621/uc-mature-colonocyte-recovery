library(ggplot2)
figdir <- "figures/stage003"
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
sourcedir <- file.path(figdir,"source_tables")
dir.create(sourcedir,recursive=TRUE,showWarnings=FALSE)
source("code/figure_theme.R")
theme003 <- function() theme_project()+theme(legend.position="none",strip.background=element_rect(fill="#F2F2F2",colour=NA),strip.text=element_text(face="bold"))
pal <- c("Clinical remission"="#2166AC","No clinical remission"="#B35806")
mean_ci <- function(x) {
 n <- sum(is.finite(x)); x <- x[is.finite(x)]
 se <- sd(x)/sqrt(n); cr <- qt(.975,n-1)
 data.frame(y=mean(x),ymin=mean(x)-cr*se,ymax=mean(x)+cr*se)
}
p <- read.delim("runs/003_cell_context/patient_metrics.tsv",check.names=FALSE,stringsAsFactors=FALSE)
p <- p[p$cell_set=="exclude_predicted_doublets" & p$scope=="all_fixed" & p$threshold==20,]
p$group <- ifelse(p$remission=="Remission","Clinical remission","No clinical remission")
make_long <- function(dat,metrics,labels){
 out <- do.call(rbind,lapply(seq_along(metrics),function(i) data.frame(
   patient=dat$patient,group=dat$group,panel=labels[i],value=dat[[metrics[i]]],stringsAsFactors=FALSE)))
 out[is.finite(out$value),]
}
a <- make_long(p,c("delta_ct_logit_epi","delta_ct_state_score"),
 c("A  CT proportion within epithelium","B  Six-gene score within CT state"))
write.table(a,file.path(sourcedir,"figure003A_patient_axes.tsv"),sep="	",quote=FALSE,row.names=FALSE)
set.seed(3003)
pa <- ggplot(a,aes(group,value,colour=group))+
 geom_hline(yintercept=0,linewidth=.35,colour="#888888",linetype=2)+
 geom_jitter(width=.12,height=0,size=2.1,alpha=.85)+
 stat_summary(fun=mean,geom="point",shape=18,size=4,colour="black")+
 stat_summary(fun.data=mean_ci,geom="errorbar",width=.14,linewidth=.65,colour="black")+
 facet_wrap(~panel,scales="free_y",nrow=1)+
 scale_colour_manual(values=pal)+
 labs(x=NULL,y="Patient-level Post - Pre change",
      title="Composition and within-state expression give different signals",
      subtitle="Adalimumab cohort; each point is one patient; diamond and bars show mean and 95% t interval")+
 theme003()+theme(axis.text.x=element_text(angle=18,hjust=1))
save_panel(pa,"figure003A_composition_vs_state",9.2,4.6)

j <- read.delim("runs/003_cell_context/joint_primary_patient_metrics.tsv",check.names=FALSE,stringsAsFactors=FALSE)
j$group <- ifelse(j$remission=="Remission","Clinical remission","No clinical remission")
b <- make_long(j,c("composition_contribution","within_state_contribution","total_epi_linear_change"),
 c("A  Composition contribution","B  Within-group contribution","C  Total epithelial change"))
write.table(b,file.path(sourcedir,"figure003B_decomposition.tsv"),sep="	",quote=FALSE,row.names=FALSE)
set.seed(3004)
pb <- ggplot(b,aes(group,value,colour=group))+
 geom_hline(yintercept=0,linewidth=.35,colour="#888888",linetype=2)+
 geom_jitter(width=.12,height=0,size=2.0,alpha=.85)+
 stat_summary(fun=mean,geom="point",shape=18,size=3.8,colour="black")+
 stat_summary(fun.data=mean_ci,geom="errorbar",width=.14,linewidth=.6,colour="black")+
 facet_wrap(~panel,scales="free_y",nrow=1)+
 scale_colour_manual(values=pal)+
 labs(x=NULL,y="Patient-level Post - Pre contribution",
      title="The epithelial score difference is carried by within-group expression",
      subtitle="Jointly evaluable patients (n=16); decomposition is algebraic and not causal")+
 theme003()+theme(axis.text.x=element_text(angle=18,hjust=1),panel.spacing.x=grid::unit(3.8,"lines"))
save_panel(pb,"figure003B_kitagawa_decomposition",12.8,4.8)

# Sensitivity source is generated directly from the branch field.
m <- read.delim("runs/003_cell_context/patient_group_models.tsv",check.names=FALSE,stringsAsFactors=FALSE)
z <- m[m$metric %in% c("delta_ct_logit_epi","delta_ct_state_score") &
       grepl("^exclude_predicted_doublets::",m$branch),]
parts <- do.call(rbind,strsplit(z$branch,"::",fixed=TRUE))
z$scope <- parts[,2]; z$threshold <- sub("threshold","",parts[,3])
z$axis <- ifelse(z$metric=="delta_ct_logit_epi","CT proportion within epithelium","Six-gene score within CT state")
z$branch_label <- paste(c(all_fixed="All fixed pairs",v31_only="v3.1 only",author_both="Author-paired subset")[z$scope],
                        paste0("threshold ",z$threshold),sep=", ")
write.table(z,file.path(sourcedir,"figure003C_sensitivity_models.tsv"),sep="	",quote=FALSE,row.names=FALSE)
z$branch_label <- factor(z$branch_label,levels=rev(unique(z$branch_label)))
pc <- ggplot(z,aes(estimate,branch_label))+
 geom_vline(xintercept=0,linetype=2,colour="#777777",linewidth=.4)+
 geom_errorbar(aes(xmin=lower_HC3,xmax=upper_HC3),width=.18,linewidth=.55,colour="#444444",orientation="y")+
 geom_point(aes(fill=scope),shape=21,size=2.5,colour="black")+
 facet_wrap(~axis,scales="free_x",nrow=1)+
 scale_fill_manual(values=c(all_fixed="#2166AC",v31_only="#4DAF4A",author_both="#984EA3"),
                    breaks=c("all_fixed","v31_only","author_both"),
                    labels=c("All fixed pairs","v3.1 only","Author-paired subset"))+
 labs(x="Remission - non-remission difference (HC3 95% CI)",y=NULL,
      title="Sensitivity branches preserve the axis contrast",
      subtitle="Thresholds and subsets were fixed before expression-matrix extraction")+
 theme003()+theme(legend.position="bottom",legend.title=element_blank(),panel.spacing.x=grid::unit(4.2,"lines"))
save_panel(pc,"figure003C_sensitivity_forest",12.2,6.0)
cat("FIGURE003_COMPLETE")
