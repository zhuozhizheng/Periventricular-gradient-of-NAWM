#Step 1 statistical analysis
rm(list=ls())
library(stringr)
library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(boot)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(robustbase)
library(multcomp)
library(survival)
library(survminer)


savepath<-"D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/NDI_neg_results"

# Str<-'NDI_gene_neg';
# datapath<-"F:/MSBio/data/out_result/NDI_gene_neg"
# Data_ref<-read.csv('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/NDI/BP/NDI_Neg_ALL_GO/Enrichment_GO/GO_AllLists.csv',header = TRUE)

#savepath<-"D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/NDI_pos_results"
# Str<-'NDI_gene_pos';
# datapath<-"F:/MSBio/data/out_result/NDI_gene_pos"
# Data_ref<-read.csv('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/NDI/BP/NDI_Pos_ALL_GO/Enrichment_GO/GO_AllLists.csv',header = TRUE)

# savepath<-"D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/ODI_neg_results"
# Str<-'ODI_gene_neg';
# datapath<-"F:/MSBio/data/out_result/ODI_gene_neg"
# Data_ref<-read.csv('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Neg_ALL_GO/Enrichment_GO/GO_AllLists.csv',header = TRUE)

# savepath<-"D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/ODI_pos_results"
# Str<-'ODI_gene_pos';
# datapath<-"F:/MSBio/data/out_result/ODI_gene_pos"
# Data_ref<-read.csv('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Pos_ALL_GO/Enrichment_GO/GO_AllLists.csv',header = TRUE)

# savepath<-"D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/WMH_neg_results"
# Str<-'WMH_gene_neg';
# datapath<-"F:/MSBio/data/out_result/WMH_gene_neg"
# Data_ref<-read.csv('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/Lesion_volume/BP/Multi_gene_list_PLS_Neg_GO/Enrichment_GO/GO_AllLists.csv',header = TRUE)

savepath<-"D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/GO_permutation/WMH_pos_results"
Str<-'WMH_gene_pos';
datapath<-"F:/MSBio/data/out_result/WMH_gene_pos"
Data_ref<-read.csv('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/Lesion_volume/BP/Multi_gene_list_PLS_Pos_GO/Enrichment_GO/GO_AllLists.csv',header = TRUE)

setwd(datapath);
files<-dir(path=datapath,pattern='Gene_list*')

Item<-c('GO Biological Processes','KEGG Pathway','Hallmark Gene Sets','Immunologic Signatures')
#Ind_item<-1

HC<-NULL
AD<-NULL
PD<-NULL
SVD<-NULL
MS<-NULL



for(i in 1:1000) #length(files)
{
print(i)
  
if (file.exists(paste0(datapath,'/',files[i],'/Enrichment_GO/GO_AllLists.csv')))
{
Data<-read.csv(paste0(datapath,'/',files[i],'/Enrichment_GO/GO_AllLists.csv'),header = TRUE)

HC<-rbind(HC,Data[Data$GeneList=='HC',]);
AD<-rbind(HC,Data[Data$GeneList=='AD',]);
PD<-rbind(HC,Data[Data$GeneList=='PD',]);
SVD<-rbind(HC,Data[Data$GeneList=='SVD',]);
MS<-rbind(HC,Data[Data$GeneList=='MS',]);
}
}
colnames(HC)<-colnames(Data);
colnames(AD)<-colnames(Data);
colnames(PD)<-colnames(Data);
colnames(SVD)<-colnames(Data);
colnames(MS)<-colnames(Data);
Data_ref$X.GeneInGOAndHitList

Data_ref_new<-Data_ref
for(i in 1:dim(Data_ref)[1])
{
  tem_go<-Data_ref[i,'GO']
  print(tem_go)
  if(Data_ref[i,'GeneList']=='HC')
  {
    match_num<-length(which(HC$GO==tem_go))
    pvalue<-match_num/1000
    Data_ref_new[i,'Occure_frequency']<-match_num
   # Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
    tem_Data<-HC[which(HC$GO==tem_go,arr.ind = TRUE),]
    match_num<-length(which(tem_Data$X.GeneInGOAndHitList>Data_ref[i,'X.GeneInGOAndHitList'],arr.ind = TRUE))
    pvalue<-(match_num+1)/1000
    Data_ref_new[i,'Rank']<-match_num+1
    Data_ref_new[i,'uncorrected_per_p']<-pvalue
    # 
    # tem_Data<-HC[which(HC$GO==tem_go,arr.ind = TRUE),]
    # match_num<-length(which(tem_Data$Z.score>Data_ref[i,'Z.score'],arr.ind = TRUE))
    # pvalue<-match_num/length(files)
    # Data_ref_new[i,'Occure_frequency']<-match_num
    # Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
  }
  

  
  if(Data_ref[i,'GeneList']=='AD')
  {
    match_num<-length(which(AD$GO==tem_go))
    pvalue<-match_num/1000
    Data_ref_new[i,'Occure_frequency']<-match_num
    #Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
    tem_Data<-AD[which(AD$GO==tem_go,arr.ind = TRUE),]
    match_num<-length(which(tem_Data$X.GeneInGOAndHitList>Data_ref[i,'X.GeneInGOAndHitList'],arr.ind = TRUE))
    pvalue<-(match_num+1)/1000
    Data_ref_new[i,'Rank']<-match_num+1
    Data_ref_new[i,'uncorrected_per_p']<-pvalue
    # 
    # tem_Data<-HC[which(HC$GO==tem_go,arr.ind = TRUE),]
    # match_num<-length(which(tem_Data$Z.score>Data_ref[i,'Z.score'],arr.ind = TRUE))
    # pvalue<-match_num/length(files)
    # Data_ref_new[i,'Occure_frequency']<-match_num
    # Data_ref_new[i,'uncorrected_per_p']<-pvalue
  }
  
  
  if(Data_ref[i,'GeneList']=='PD')
  {
    match_num<-length(which(PD$GO==tem_go))
    pvalue<-match_num/1000
    Data_ref_new[i,'Occure_frequency']<-match_num
   # Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
    tem_Data<-PD[which(PD$GO==tem_go,arr.ind = TRUE),]
    match_num<-length(which(tem_Data$X.GeneInGOAndHitList>Data_ref[i,'X.GeneInGOAndHitList'],arr.ind = TRUE))
    pvalue<-(match_num+1)/1000
    Data_ref_new[i,'Rank']<-match_num+1
    Data_ref_new[i,'uncorrected_per_p']<-pvalue
    # 
    # tem_Data<-HC[which(HC$GO==tem_go,arr.ind = TRUE),]
    # match_num<-length(which(tem_Data$Z.score>Data_ref[i,'Z.score'],arr.ind = TRUE))
    # pvalue<-match_num/length(files)
    # Data_ref_new[i,'Occure_frequency']<-match_num
    # Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
  }
  
  
  if(Data_ref[i,'GeneList']=='SVD')
  {
    match_num<-length(which(SVD$GO==tem_go))
    pvalue<-match_num/1000
    Data_ref_new[i,'Occure_frequency']<-match_num
    #Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
    tem_Data<-SVD[which(SVD$GO==tem_go,arr.ind = TRUE),]
    match_num<-length(which(tem_Data$X.GeneInGOAndHitList>Data_ref[i,'X.GeneInGOAndHitList'],arr.ind = TRUE))
    pvalue<-(match_num+1)/1000
    Data_ref_new[i,'Rank']<-match_num+1
    Data_ref_new[i,'uncorrected_per_p']<-pvalue
    # 
    # tem_Data<-HC[which(HC$GO==tem_go,arr.ind = TRUE),]
    # match_num<-length(which(tem_Data$Z.score>Data_ref[i,'Z.score'],arr.ind = TRUE))
    # pvalue<-match_num/length(files)
    # Data_ref_new[i,'Occure_frequency']<-match_num
    # Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
  }
  
  
  if(Data_ref[i,'GeneList']=='MS')
  {
    match_num<-length(which(MS$GO==tem_go))
    pvalue<-match_num/1000
    Data_ref_new[i,'Occure_frequency']<-match_num
    #Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
    tem_Data<-MS[which(MS$GO==tem_go,arr.ind = TRUE),]
    match_num<-length(which(tem_Data$X.GeneInGOAndHitList>Data_ref[i,'X.GeneInGOAndHitList'],arr.ind = TRUE))
    pvalue<-(match_num+1)/1000
    Data_ref_new[i,'Rank']<-match_num+1
    Data_ref_new[i,'uncorrected_per_p']<-pvalue
    # 
    # tem_Data<-HC[which(HC$GO==tem_go,arr.ind = TRUE),]
    # match_num<-length(which(tem_Data$Z.score>Data_ref[i,'Z.score'],arr.ind = TRUE))
    # pvalue<-match_num/length(files)
    # Data_ref_new[i,'Occure_frequency']<-match_num
    # Data_ref_new[i,'uncorrected_per_p']<-pvalue
    
  }
  
}


Data_ref_new[,'pFDR_per']<-p.adjust(Data_ref_new$uncorrected_per_p,'fdr')


Item<-c('GO Biological Processes','KEGG Pathway','Hallmark Gene Sets','Immunologic Signatures')

Data_ref_new[Data_ref_new$GeneList=='GO Biological Processes','pFDR_per']<-p.adjust(Data_ref_new[Data_ref_new$GeneList=='GO Biological Processes','uncorrected_per_p'],'fdr')
Data_ref_new[Data_ref_new$GeneList=='KEGG Pathway','pFDR_per']<-p.adjust(Data_ref_new[Data_ref_new$GeneList=='KEGG Pathway','uncorrected_per_p'],'fdr')
Data_ref_new[Data_ref_new$GeneList=='Hallmark Gene Sets','pFDR_per']<-p.adjust(Data_ref_new[Data_ref_new$GeneList=='Hallmark Gene Sets','uncorrected_per_p'],'fdr')
Data_ref_new[Data_ref_new$GeneList=='Immunologic Signatures','pFDR_per']<-p.adjust(Data_ref_new[Data_ref_new$GeneList=='Immunologic Signatures','uncorrected_per_p'],'fdr')

setwd(savepath)
write.csv(Data_ref_new,paste0(Str,'_Metascape_output.csv'));

# #arrnage the data for heatmap 
# gene_ref <- arrange(gene_ref,gene_ref$pFDR_per)
# 
# 
# setwd(savepath)