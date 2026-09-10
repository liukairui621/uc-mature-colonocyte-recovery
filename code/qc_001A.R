options(stringsAsFactors=FALSE, width=150)
Sys.setenv(OMP_NUM_THREADS=4,OPENBLAS_NUM_THREADS=4)
suppressPackageStartupMessages({library(limma);library(data.table);library(jsonlite);library(ggplot2)})
root <- "/root/projects/UC_Treatment_Recovery"
out <- file.path(root,"runs/001A_qc")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20260910)
input <- file.path(root,"inputs/discovery/GSE92415_series_matrix.relay.txt.gz")
read_geo <- function(path) {
 con<-gzfile(path,"rt");on.exit(close(con))
 repeat{line<-readLines(con,n=1,warn=FALSE);if(!length(line))stop("No matrix start");if(line=="!series_matrix_table_begin")break}
 d<-read.delim(con,check.names=FALSE,quote='"',comment.char="!",row.names=1)
 as.matrix(d)
}
read_annot <- function(path) {
 con<-gzfile(path,"rt");on.exit(close(con))
 repeat{line<-readLines(con,n=1,warn=FALSE);if(!length(line))stop("No annotation table");if(startsWith(line,"ID\t"))break}
 header<-strsplit(line,"\t",fixed=TRUE)[[1]]
 d<-read.delim(con,header=FALSE,col.names=header,check.names=FALSE,quote="",comment.char="!",fill=TRUE)
 d
}
x<-read_geo(input)
storage.mode(x)<-"double"
m<-fread(file.path(root,"runs/001A_metadata/sample_registry_v1.tsv"),na.strings="")
m<-as.data.frame(m[m$series=="GSE92415",]);m<-m[match(colnames(x),m$gsm),]
stopifnot(nrow(x)==54715,ncol(x)==183,identical(m$gsm,colnames(x)),!anyDuplicated(rownames(x)),all(is.finite(x)))
q<-quantile(x,probs=c(0,.01,.25,.5,.75,.99,1))
# Source documents RMA normalization; the numerical range must be compatible with log2.
if(q[7]>40 || q[1]< -20)stop("Scale inconsistent with expected RMA log2; inspect before transformation")
a<-read_annot(file.path(root,"inputs/annotation/GPL13158.annot.gz"))
stopifnot("Gene symbol"%in%names(a))
ai<-match(rownames(x),a$ID);sy<-a[["Gene symbol"]][ai]
keep<-!is.na(sy)&nzchar(sy)&!grepl("///|//|;|\\|",sy)&!grepl("^AFFX",rownames(x))
mapping<-data.frame(probe=rownames(x),gene_symbol=sy,retained=keep)
fwrite(mapping,file.path(out,"probe_gene_mapping.tsv"),sep="\t",na="NA")
g<-avereps(x[keep,,drop=FALSE],ID=sy[keep])
g<-g[order(rownames(g)),,drop=FALSE]
stopifnot(!anyDuplicated(rownames(g)),all(is.finite(g)))
qs<-t(apply(x,2,quantile,probs=c(0,.25,.5,.75,1)))
qc<-data.frame(gsm=colnames(x),missing=colSums(!is.finite(x)),minimum=qs[,1],q25=qs[,2],median=qs[,3],q75=qs[,4],maximum=qs[,5],sd=apply(x,2,sd),m[,c("subject","disease","treatment","time","response")])
qc$group<-ifelse(qc$disease=="Healthy","Healthy",paste(qc$treatment,qc$time,sep="_"))
vars<-apply(g,1,var);top<-order(vars,decreasing=TRUE)[seq_len(min(5000,nrow(g)))]
pc<-prcomp(t(g[top,,drop=FALSE]),center=TRUE,scale.=FALSE,rank.=10)
pv<-pc$sdev^2/sum(pc$sdev^2)
pcs<-data.frame(qc,pc$x);fwrite(pcs,file.path(out,"sample_qc_and_pca.tsv"),sep="\t",na="NA")
fwrite(data.frame(component=seq_along(pv),variance_fraction=pv),file.path(out,"pca_variance.tsv"),sep="\t")
# No biological outlier is automatically removed. Exact duplicates and nonfinite data are checked.
dups<-which(duplicated(as.data.frame(t(x))))
cors<-cor(g)
diag(cors)<-NA_real_
fwrite(data.frame(gsm=colnames(g),max_other_correlation=apply(cors,2,max,na.rm=TRUE),median_other_correlation=apply(cors,2,median,na.rm=TRUE)),file.path(out,"sample_correlations.tsv"),sep="\t")
pal<-c(Healthy="#777777",golimumab_Week_0="#B8DCE8",golimumab_Week_6="#167C9B",Placebo_Week_0="#F3D1B2",Placebo_Week_6="#B45E2B")
qc$plot_group<-gsub(" ","_",qc$group);pcs$plot_group<-gsub(" ","_",pcs$group)
theme_uc<-function(){theme_classic(base_size=11,base_family="sans")+theme(legend.position="bottom",plot.title=element_text(face="bold"),axis.line=element_line(linewidth=.4),axis.ticks=element_line(linewidth=.4))}
save_plot<-function(p,name,w=8,h=5.6){ggsave(file.path(out,paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf);ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=160)}
p<-ggplot(pcs,aes(PC1,PC2,color=plot_group))+geom_point(size=2,alpha=.8)+scale_color_manual(values=pal,name=NULL)+labs(title="Discovery cohort: expression PCA",subtitle="Top 5,000 variable genes; no automatic sample exclusion",x=sprintf("PC1 (%.1f%%)",100*pv[1]),y=sprintf("PC2 (%.1f%%)",100*pv[2]))+theme_uc()
save_plot(p,"QC01_PCA")
p<-ggplot(qc,aes(x=plot_group,y=median,color=plot_group))+geom_boxplot(outlier.shape=NA,width=.6)+geom_jitter(width=.12,height=0,size=1.5,alpha=.65)+scale_color_manual(values=pal,guide="none")+labs(title="Sample median RMA expression",x=NULL,y="Median log2 expression")+theme_uc()+theme(axis.text.x=element_text(angle=20,hjust=1))
save_plot(p,"QC02_Medians")
pairs<-as.data.frame(fread(file.path(root,"runs/001A_metadata/discovery_pairs_v1.tsv")))
stopifnot(all(pairs$baseline_gsm%in%colnames(g)),all(pairs$post_gsm%in%colnames(g)),!anyDuplicated(pairs$subject))
delta<-g[,pairs$post_gsm,drop=FALSE]-g[,pairs$baseline_gsm,drop=FALSE];colnames(delta)<-pairs$subject
saveRDS(list(probe_expression=x,gene_expression=g,metadata=m,pairs=pairs,delta=delta,annotation=mapping,quantiles=q),file.path(out,"discovery_preprocessed_v1.rds"),compress="gzip")
summary<-list(run_id="001A_qc",stage="technical QC only; no candidate frozen",status="COMPLETE",input=input,probes=nrow(x),arrays=ncol(x),mapped_single_gene_probes=sum(keep),genes=nrow(g),quantiles=as.list(q),nonfinite_values=sum(!is.finite(x)),constant_probes=sum(apply(x,1,var)==0),exact_duplicate_columns=length(dups),duplicate_column_names=colnames(x)[dups],pairs=nrow(pairs),pairs_by_arm=as.list(table(pairs$arm)),normalization="Source ArrayStudio v9 RMA retained; no second log2 or renormalization",probe_collapse="Arithmetic mean of log2 probe values per uniquely mapped gene; ambiguous symbols and AFFX controls excluded",sample_exclusions="No QC exclusions applied; primary analysis restricted to verified paired subjects",healthy_comparison="Descriptive within-cohort reference; separate recruitment from trial",validation_expression_accessed=FALSE)
write_json(summary,file.path(out,"qc_summary.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat(toJSON(summary,pretty=TRUE,auto_unbox=TRUE),"\n")
