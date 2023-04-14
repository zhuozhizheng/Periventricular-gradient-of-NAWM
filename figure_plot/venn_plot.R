rm(list=ls())
library(xlsx)
library(VennDiagram)
path<-'D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI'
savepath<-'D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP'
setwd(path)
list1<-list(1,5)


HC<-read.csv('HC_Genename_column1.csv',header=FALSE)
AD<-read.csv('AD_Genename_column1.csv',header=FALSE)
PD<-read.csv('PD_Genename_column1.csv',header=FALSE)
SVD<-read.csv('SVD_Genename_column1.csv',header=FALSE)
MS<-read.csv('MS_Genename_column1.csv',header=FALSE)


# inter1<-intersect(HC$V1,AD$V1)
# inter2<-intersect(inter1,PD$V1)
# inter3<-intersect(inter2,SVD$V1)
# inter4<-intersect(inter3,MS$V1)
# 
# inter4<-data.frame(inter4)


setwd(savepath)

#for PLS-
 # HC1<-HC[HC$V2<0,]
 # AD1<-AD[AD$V2<0,]
 # PD1<-PD[PD$V2<0,]
 # SVD1<-SVD[SVD$V2<0,]
 # MS1<-MS[MS$V2<0,]
 # 
 # inter1<-intersect(HC1$V1,AD1$V1)
 # inter2<-intersect(inter1,PD1$V1)
 # inter3<-intersect(inter2,SVD1$V1)
 # inter4<-intersect(inter3,MS1$V1)
 # 
 # inter4<-data.frame(inter4)
 # write.csv(inter4,'Multi_gene_list_PLS_Neg_All_disease_overlap.csv')

 HC1<-HC[HC$V2>0,]
 AD1<-AD[AD$V2>0,]
 PD1<-PD[PD$V2>0,]
 SVD1<-SVD[SVD$V2>0,]
 MS1<-MS[MS$V2>0,]

 inter1<-intersect(HC1$V1,AD1$V1)
 inter2<-intersect(inter1,PD1$V1)
 inter3<-intersect(inter2,SVD1$V1)
 inter4<-intersect(inter3,MS1$V1)

 inter4<-data.frame(inter4)
write.csv(inter4,'Multi_gene_list_PLS_Pos_All_disease_overlap.csv')

colnames(HC1)[1]<-c('z1')
colnames(AD1)[1]<-c('z2')
colnames(PD1)[1]<-c('z3')
colnames(SVD1)[1]<-c('z4')
colnames(MS1)[1]<-c('z5')

#list_data<-data.frame(HC1$V1,AD1$V1,PD1$V1,SVD1$V1,MS1$V1)


p <- venn.diagram(
  list(HC=HC1$z1,AD=AD1$z2,PD=PD1$z3,SVD=SVD1$z4,MS=MS1$z5),
  filename = "out5venn_pos.tiff",
  lty = "dotted",
  lwd = 2,
  col = "black",  #"transparent",
  #fill = c("dodgerblue", "goldenrod1", "darkred1", "seagreen3", "orchid3"),
  alpha = 1,
  #cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "green3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.1,
  cex = 1.5
)

#inter<-get.venn.partitions(list_data)

print(p)

#dev.off()