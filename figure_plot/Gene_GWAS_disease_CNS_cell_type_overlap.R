#library(clusterProfiler)
#library(DOSE)

#library(org.Hs.eg.db)
#library(topGO)
#library(pathview)
library(ggplot2)
#library(readxl)
#library(enrichplot)
library(Cairo)
rm(list=ls())


setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/Lesion_volume')

 # gene1<-read.csv('Cell_Neg.csv',header=TRUE)
 # str='Cell_Neg'
gene1<-read.csv('Cell_Pos.csv',header=TRUE)
str='Cell_Pos'
# 
 # gene1<-read.csv('Neg_disease_result.csv',header=TRUE)
 # str='Neg_disease_result'
# gene1<-read.csv('Pos_disease_result.csv',header=TRUE)
# str='Pos_disease_result'

 

colnames(gene1)[1]<-'Description'
gene1$GeneList <-factor(gene1$GeneList,levels=c('HC','AD','PD','SVD','MS'))



gene1$Logqvalue<-log10(gene1$pFDR)
 
gene1$Logqvalue<-(-gene1$Logqvalue)
#gene1$Logqvalue<-gene1$pFDR
gene1$Ratio<-gene1$InGO#/gene1$Hitlist

lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 16, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))

#pdf(file=paste0(str,"_Popo_figure.pdf"),width=15,height=10)

p<-ggplot2::ggplot(gene1,aes(x=GeneList,y=reorder(Description,Ratio),size=Ratio,colour=Logqvalue))+
  
  geom_point(shape=16)+
  labs(X=" ",y=" ")+
  scale_colour_continuous(name="-log(pFDR)", low="blue", high="red")+
  scale_radius(range = c(0,10),name="GeneNumber")+

  #guides(color=guide_colorbar(order=1),size=guide_legend(order=2))+
  #guides(color=guide_colorbar(order=1),size=guide_legend(order=2))+
  theme(text = element_text(family = "Arial",size=16,colour='black'),
        axis.text = element_text(colour='black'),
        axis.text.x = element_text(colour='black'),
        axis.text.y = element_text(colour='black'),
        axis.title.x = element_text(colour='black'),
        axis.title.y = element_text(colour='black'),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))

  #theme_bw()
  #theme_bw(base_size = 20, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt)

#print(p1)
#dev.off()

ggsave(paste0(str,"_Rplot01.pdf"),p,device = cairo_pdf,width=15,height=12,units='cm')