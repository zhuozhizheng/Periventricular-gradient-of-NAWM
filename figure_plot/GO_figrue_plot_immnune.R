library(clusterProfiler)
library(DOSE)

library(org.Hs.eg.db)
library(topGO)
library(pathview)
library(ggplot2)
library(readxl)
library(enrichplot)
library(Cairo)
rm(list=ls())
# 
# Str<-'NDI_neg_results';
# savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/NDI_neg_results')
# gene3<-read.csv('NDI_gene_neg_Metascape_output.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/NDI/BP/NDI_Neg_ALL_Immune/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGOTop100.csv',header=TRUE)

# Str<-'NDI_pos_results';
# savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/NDI_pos_results')
# gene3<-read.csv('NDI_gene_pos_Metascape_output.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/NDI/BP/NDI_Pos_ALL_Immune/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGOTop100.csv',header=TRUE)


# Str<-'ODI_neg_results';
# savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/ODI_neg_results')
# gene3<-read.csv('ODI_gene_neg_Metascape_output.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Neg_ALL_GO_Immune/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGOTop100.csv',header=TRUE)

Str<-'ODI_pos_results';
savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/ODI_pos_results')
gene3<-read.csv('ODI_gene_pos_Metascape_output.csv',header=TRUE)
setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Pos_ALL_GO_Immune/Enrichment_heatmap')
gene_ref<-read.csv('HeatmapSelectedGOTop100.csv',header=TRUE)


# Str<-'WMH_neg_results';
# savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/WMH_neg_results')
# gene3<-read.csv('WMH_gene_neg_Metascape_output.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/Lesion_volume/BP/Multi_gene_list_PLS_Neg_immune/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGOTop100.csv',header=TRUE)

# Str<-'WMH_pos_results';
# savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/WMH_pos_results')
# gene3<-read.csv('WMH_gene_pos_Metascape_output.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/Lesion_volume/BP/Multi_gene_list_PLS_Pos_immune/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGOTop100.csv',header=TRUE)

gene3<-gene3 %>% filter(Category=='Immunologic Signatures')
gene3<-gene3 %>% filter(pFDR_per<0.05)
gene1<-gene3
#gene1[gene1$LogP>-3,'LogP']<-NA

# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Pos_ALL_GO/Enrichment_GO')
# gene1<-read.csv('GO_AllLists.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Pos_ALL_GO/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)


# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan/Results/NDI/BP/PLS-Endo/Enrichment_GO')
# gene1<-read.csv('GO_AllLists.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan/Results/NDI/BP/PLS-Endo/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)


#setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan/Results/NDI/BP/Multi_gene_list_PLS_Pos_immune/Enrichment_GO')



#dotplot(gene1,showCategory=10,includeAll=TRUE)
#gene1$Log.q.value.<-(-gene1$Log.q.value.)
gene1$LogP<-(-gene1$LogP)
gene1$Ratio<-gene1$X.GeneInGOAndHitList#/gene1$X.GeneInHitList
# gene1<-gene1 %>% filter(Log.q.value.>0)

#gene1<-gene2

#index_sel<-intersect(gene1$Description,gene2$Description)
 
 index<-NULL
 for (i in 1:nrow(gene1))
 {
 
   if (!is.na(match(gene1$Description[i],gene_ref$Description)))
   {index<-c(index,i)}
 }
 
 gene1<-gene1[index,]

  
# for (i in 1:nrow(gene1))
# {
#   if (!length(setdiff(gene1$Category[i],'GO Biological Processes')))
#   {gene1[i,'Description']<-paste0('GO:',gene1$Description[i])}
#   else if (!length(setdiff(gene1$Category[i],'KEGG Pathway')))
#   {gene1[i,'Description']<-paste0('KEGG:',gene1$Description[i])}
#   else if (!length(setdiff(gene1$Category[i],'Hallmark Gene Sets')))
#   {gene1[i,'Description']<-paste0('H:',gene1$Description[i])}
#   
# }

# for (i in 1:nrow(gene1))
# {
#   if (!length(setdiff(gene1$GeneList[i],'A_HC')))
#   {gene1[i,'GeneList']<-paste0('HC')}
#   else if (!length(setdiff(gene1$GeneList[i],'B_AD')))
#   {gene1[i,'GeneList']<-paste0('AD')}
#   else if (!length(setdiff(gene1$GeneList[i],'C_PD')))
#   {gene1[i,'GeneList']<-paste0('PD')}
#   else if (!length(setdiff(gene1$GeneList[i],'D_SVD')))
#   {gene1[i,'GeneList']<-paste0('SVD')}
#   else if (!length(setdiff(gene1$GeneList[i],'F_MS')))
#   {gene1[i,'GeneList']<-paste0('MS')}
#   
# }

#gene1$GeneList1<-NULL
# gene1$GeneList[A_HC]<-'HC'
# gene1$GeneList[B_AD]<-'AD'
# gene1$GeneList[A_PD]<-'PD'
# gene1$GeneList[A_SVD]<-'SVD'
# gene1$GeneList[A_MS]<-'MS'
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 16, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))


#pdf(file="Popo_figure.pdf")

gene1$GeneList <-factor(gene1$GeneList,levels=c('HC','AD','PD','SVD','MS'))


#gene1[,'log_pFDR_per']<--log(gene1$pFDR_per)


#Sum_value<-aggregate(gene1$Log.q.value.,by=list(type=gene1$Description),sum)

Sum_value<-aggregate(gene1$LogP,by=list(type=gene1$Description),sum)

y_label<-reorder(Sum_value$type,Sum_value$x)

#p<-ggplot(gene1,aes(x=GeneList,y=reorder(Description,Log.q.value.),size=Ratio,colour=Log.q.value.))+
p<-ggplot(gene1,aes(x=GeneList,y=reorder(Description,LogP),size=Ratio,colour=LogP))+
#p<-ggplot(gene1,aes(x=GeneList,y=reorder(Description,log_pFDR_per),size=Ratio,colour=log_pFDR_per))+
#p<-ggplot(gene1,aes(x=GeneList,y=Description,size=Ratio,colour=Log.q.value.))+
  geom_point(shape=16)+
  labs(X=" ",y=" ")+
  scale_colour_continuous(name="-log(p)", low="blue", high="red")+
  scale_radius(range = c(1,6),name="GeneNumber")+
  
  # theme(axis.text   = element_text(size = rel(0.8)), 
  #       strip.text  = element_text(size = rel(0.8)),
  #       legend.text = element_text(size = rel(0.8)),
  #       plot.title  = element_text(size = rel(1.2)),
  #       panel.grid.minor = element_line(size = rel(0.5)))+
  # guides(color=guide_colorbar(order=1),size=guide_legend(order=2))+
  theme(text = element_text(family = "Arial",
                            size=14,colour='black'),
        axis.text = element_text(colour='black'),
        axis.text.x = element_text(colour='black'),
        axis.text.y = element_text(colour='black'),
        axis.title.x = element_text(colour='black'),
        axis.title.y = element_text(colour='black'),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12)) +
   scale_y_discrete(limits = levels(y_label)[(length(y_label)-20):length(y_label)])


 # theme_bw(base_size = 20, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt)
#print(p)
#dev.off()
setwd(savepath);
ggsave(paste0(Str,"_Rplot01_Immune.pdf"),p,device = cairo_pdf,width=30,height=12,units='cm')