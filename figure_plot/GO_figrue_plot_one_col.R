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



# Str<-'NDI_neg_results';
# savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# # 
# 
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/NDI/BP/Multi_gene_list_PLS_Neg_All_disease_overlap/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/NDI_neg_results')
# gene3<-read.csv('NDI_gene_neg_Metascape_output.csv',header=TRUE)



Str<-'NDI_pos_results';
savepath=paste0('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/',Str);
# 

setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/NDI/BP/Multi_gene_list_PLS_pos_All_disease_overlap/Enrichment_heatmap')
gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)
setwd('D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/NDI_pos_results')
gene3<-read.csv('NDI_gene_pos_Metascape_output.csv',header=TRUE)


gene3<-gene3 %>% filter(Category!='Immunologic Signatures')
gene3<-gene3 %>% filter(pFDR_per>=0.05)

#gene1<-gene3

index_del<-NULL

for (i in 1:dim(gene_ref)[1])
  {
  
  
  if (length(which(gene3$Description==gene_ref$Description[i],arr.ind = TRUE))>0)
      {
        
        index_del<-c(index_del,which(gene3$Description==gene_ref$Description[i],arr.ind = TRUE))
        
  }
}


gene_ref<-gene_ref[-index_del,]
# #dotplot(gene1,showCategory=10,includeAll=TRUE)
# gene1$Log.q.value.<-(-gene1$Log.q.value.)
# gene1$Ratio<-gene1$X.GeneInGOAndHitList#/gene1$X.GeneInHitList
# gene2<-gene1 %>% filter(Log.q.value.>0)
# 
# 
# 
# index_sel<-intersect(gene1$Description,gene2$Description)
# 
# index<-NULL
# for (i in 1:nrow(gene2))
# {
#   
#   if (!is.na(match(gene2$Description[i],gene_ref$Description)))
#   {index<-c(index,i)}
# }
# 
# gene1<-gene2[index,] 
# 
#   
for (i in 1:nrow(gene_ref))
{
  if (!length(setdiff(gene_ref$Category[i],'GO Biological Processes')))
  {gene_ref[i,'Description']<-paste0('GO:',gene_ref$Description[i])}
  else if (!length(setdiff(gene_ref$Category[i],'KEGG Pathway')))
  {gene_ref[i,'Description']<-paste0('KEGG:',gene_ref$Description[i])}
  else if (!length(setdiff(gene_ref$Category[i],'Hallmark Gene Sets')))
  {gene_ref[i,'Description']<-paste0('H:',gene_ref$Description[i])}

}
# 
# # for (i in 1:nrow(gene1))
# # {
# #   if (!length(setdiff(gene1$GeneList[i],'A_HC')))
# #   {gene1[i,'GeneList']<-paste0('HC')}
# #   else if (!length(setdiff(gene1$GeneList[i],'B_AD')))
# #   {gene1[i,'GeneList']<-paste0('AD')}
# #   else if (!length(setdiff(gene1$GeneList[i],'C_PD')))
# #   {gene1[i,'GeneList']<-paste0('PD')}
# #   else if (!length(setdiff(gene1$GeneList[i],'D_SVD')))
# #   {gene1[i,'GeneList']<-paste0('SVD')}
# #   else if (!length(setdiff(gene1$GeneList[i],'F_MS')))
# #   {gene1[i,'GeneList']<-paste0('MS')}
# #   
# # }
# 
# #gene1$GeneList1<-NULL
# # gene1$GeneList[A_HC]<-'HC'
# # gene1$GeneList[B_AD]<-'AD'
# # gene1$GeneList[A_PD]<-'PD'
# # gene1$GeneList[A_SVD]<-'SVD'
# # gene1$GeneList[A_MS]<-'MS'
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 16, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))
# 


# #按照PValue从低到高排序[升序]
# gene_ref <- arrange(gene_ref,gene_ref$LogP)
# #Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
# #将Pathway列转化为因子型
# gene_ref$Term <- factor(gene_ref$Term,levels = rev(gene_ref$Term))
# 
# dev.new()
# mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), #坐标轴标题
#                  axis.text=element_text(face="bold", size=14,colour = 'black'), #坐标轴标签
#                  axis.line = element_line(size=0.5, colour = 'black'), #轴线
#                  panel.background = element_rect(color='black'), #绘图区边框
#                  legend.key = element_blank() #关闭图例边框
# )
# 
# 
# #绘制KEGG气泡图
# p <- ggplot(gene_ref,aes(x=InTerm_InList,y=description,colour=-1*LogP,size=GeneNumber))+
#   geom_point()+
#   scale_size(range=c(2, 8))+
#   scale_colour_gradient(low = "blue",high = "red")+
#   theme_bw()+
#   ylab("KEGG Pathway Terms")+
#   xlab("Gene Ratio")+
#   labs(color=expression(-log[10](PValue)))+
#   theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
#   theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
#   theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
# plot <- p+mytheme
# plot
# 
# dev.off()
# 
# #保存图片
# ggsave(plot,filename = "KEGG.pdf",width = 10,height = 6,dpi=300)
# ggsave(plot,filename = "KEGG.png",width = 10,height = 6,dpi=300)
# 
# 
# 
# 
# ggsave(plot,filename = "KEGG.pdf",width = 10,height = 6,dpi=300)
# ggsave(plot,filename = "KEGG.png",width = 10,height = 6,dpi=300)
# 

#pdf(file="Popo_figure.pdf")
y_label=reorder(gene_ref$Description,-gene_ref$LogP)
#y_label=reorder(gene_ref$Description,gene_ref$GeneNumber/gene_ref$number)

p<-ggplot(gene_ref,aes(x=GeneNumber/number,y=reorder(Description,-LogP),size=GeneNumber,colour=-LogP))+
#p<-ggplot(gene1,aes(x=GeneList,y=Description,size=Ratio,colour=Log.q.value.))+
  geom_point(shape=16)+
  labs(X=" ",y=" ")+
  scale_colour_continuous(name="-log(p)", low="blue", high="red")+
  scale_radius(range = c(3,10),name="GeneNumber")+
  # theme(axis.text   = element_text(size = rel(0.8)), 
  #       strip.text  = element_text(size = rel(0.8)),
  #       legend.text = element_text(size = rel(0.8)),
  #       plot.title  = element_text(size = rel(1.2)),
  #       panel.grid.minor = element_line(size = rel(0.5)))+
  # guides(color=guide_colorbar(order=1),size=guide_legend(order=2))+
  theme(text = element_text(family = "Arial",
                            size=20,colour='black'),
        axis.text = element_text(colour='black'),
        axis.text.x = element_text(colour='black'),
        axis.text.y = element_text(colour='black'),
        axis.title.x = element_text(colour='black'),
        axis.title.y = element_text(colour='black')) +
  scale_y_discrete(limits = levels(y_label))
  
  #scale_y_discrete(limits = levels(y_label)[(length(y_label)-20):length(y_label)])

 # theme_bw(base_size = 20, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt)
print(p)
#dev.off()
setwd(savepath);
ggsave(paste0(Str,"_Rplot01_single_col.pdf"),p,device = cairo_pdf,width=30,height=12,units='cm')