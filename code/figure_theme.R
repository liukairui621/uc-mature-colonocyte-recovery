library(ggplot2)
project_palette<-c("golimumab"="#2166AC","Placebo"="#777777","Yes"="#2166AC","No"="#B35806")
theme_project<-function(){theme_classic(base_size=11,base_family="sans")+theme(plot.title=element_text(face="bold",size=12),plot.subtitle=element_text(size=10),axis.title=element_text(size=11),axis.text=element_text(size=10),legend.title=element_text(size=10),legend.text=element_text(size=10),strip.text=element_text(size=10),legend.position="bottom",panel.spacing.x=grid::unit(3.0,"lines"),plot.margin=margin(12,18,12,12))}
save_panel<-function(p,name,width,height){
 ggsave(file.path(figdir,paste0(name,".pdf")),p,width=width,height=height,device=cairo_pdf)
 ggsave(file.path(figdir,paste0(name,".png")),p,width=width,height=height,dpi=300)
}
