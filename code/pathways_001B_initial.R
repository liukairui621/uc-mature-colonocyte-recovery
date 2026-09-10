options(stringsAsFactors=FALSE,width=150)
suppressPackageStartupMessages({library(limma);library(data.table);library(jsonlite);library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db);library(ggplot2)})
root<-"/root/projects/UC_Treatment_Recovery";set.seed(20260910)
obj<-readRDS(file.path(root,"runs/001A_qc/discovery_preprocessed_v1.rds"))
y<-obj$delta;genes<-rownames(y);out<-file.path(root,"runs/001B_discovery")
gmt<-strsplit(readLines(file.path(root,"inputs/annotation/h.all.v2025.1.Hs.symbols.gmt")),"\t",fixed=TRUE)
hallmark<-setNames(lapply(gmt,function(x)unique(x[-c(1,2)])),sapply(gmt,"[",1))
# GO membership is version pinned to installed OrgDb and GO.db (both 3.18.0).
valid<-intersect(genes,keys(org.Hs.eg.db,keytype="SYMBOL"))
go<-AnnotationDbi::select(org.Hs.eg.db,keys=valid,columns=c("GOALL","ONTOLOGYALL"),keytype="SYMBOL")
go<-unique(go[!is.na(go$GOALL)&go$ONTOLOGYALL=="BP",c("SYMBOL","GOALL")])
go_sets<-split(go$SYMBOL,go$GOALL);go_sets<-lapply(go_sets,unique)
term<-AnnotationDbi::select(GO.db,keys=names(go_sets),columns="TERM",keytype="GOID")
term_map<-setNames(term$TERM,term$GOID)
sets<-list(Hallmark=hallmark,GO_BP=go_sets)
saveRDS(sets,file.path(root,"inputs/annotation/gene_sets_001B.rds"))
capture.output(AnnotationDbi::metadata(org.Hs.eg.db),file=file.path(root,"inputs/annotation/OrgDb_metadata.txt"))
capture.output(AnnotationDbi::metadata(GO.db),file=file.path(root,"inputs/annotation/GOdb_metadata.txt"))
summary<-list()
for(branch in c("primary_paired","barcode_sensitivity")){
 bdir<-file.path(out,branch);des<-as.data.frame(fread(file.path(bdir,"main_design.tsv")));X<-as.matrix(des[,-1,drop=FALSE])
 stopifnot(identical(des$subject,colnames(y)))
 C<-matrix(0,nrow=ncol(X),ncol=3,dimnames=list(colnames(X),c("GLM_change","Placebo_change","GLM_minus_Placebo")))
 C["GLM","GLM_change"]<-1;C["Placebo","Placebo_change"]<-1;C[c("GLM","Placebo"),"GLM_minus_Placebo"]<-c(1,-1)
 z<-as.data.frame(fread(file.path(bdir,"response_design.tsv")));Z<-as.matrix(z[,-1,drop=FALSE])
 cr<-setNames(rep(0,ncol(Z)),colnames(Z));cr["G_Yes"]<-1;cr["G_No"]<- -1
 for(lib in names(sets)){
  idx<-lapply(sets[[lib]],function(g)which(genes%in%g))
  idx<-idx[lengths(idx)>=15 & lengths(idx)<=500]
  cov<-data.frame(set=names(idx),original_size=lengths(sets[[lib]][names(idx)]),mapped_size=lengths(idx))
  cov$coverage<-cov$mapped_size/cov$original_size
  fwrite(cov,file.path(bdir,paste0(lib,"_coverage.tsv")),sep="\t")
  ans<-list()
  for(co in c(colnames(C),"GLM_response_difference")){
   xx<-if(co=="GLM_response_difference")Z else X
   cc<-if(co=="GLM_response_difference")cr else C[,co]
   # Competitive gene-set test with correlation adjustment, using all tested genes.
   rr<-camera(y,index=idx,design=xx,contrast=cc,inter.gene.cor=.01,sort=TRUE,use.ranks=FALSE)
   rr$set<-rownames(rr);rr$term<-if(lib=="GO_BP")unname(term_map[rr$set]) else sub("^HALLMARK_","",rr$set)
   rr$mean_log2_change<-vapply(idx[rr$set],function(ii){
      tt<-fread(file.path(bdir,paste0(co,".tsv")))
      mean(tt$logFC[match(genes[ii],tt$gene)])
   },numeric(1))
   # Per contrast/library FDR is the primary screening output; global FDR below.
   rr$contrast<-co;rr$library<-lib;rr$branch<-branch
   fwrite(rr,file.path(bdir,paste0("camera_",lib,"_",co,".tsv")),sep="\t")
   ans[[co]]<-rr
  }
  allres<-rbindlist(ans);allres$FDR_across_four_contrasts<-p.adjust(allres$PValue,"BH")
  fwrite(allres,file.path(bdir,paste0("camera_",lib,"_all.tsv")),sep="\t")
  summary[[paste(branch,lib,sep="_")]]<-as.list(table(allres$contrast[allres$FDR<.05]))
 }
}
write_json(list(status="COMPLETE",stage="exploratory program screen",method="limma camera on patient-paired changes; competitive test, fixed inter-gene correlation 0.01; ranks=FALSE",libraries=c("MSigDB Hallmark 2025.1.Hs","GO BP from org.Hs.eg.db/GO.db 3.18.0"),sizes=c(min=15,max=500),counts=summary,interpretation="Relative coordinated transcriptional changes; not direct functional activity or causal drug specificity",candidate_frozen=FALSE),file.path(out,"pathway_summary.json"),auto_unbox=TRUE,pretty=TRUE)
cat(toJSON(summary,auto_unbox=TRUE,pretty=TRUE),"\n")
