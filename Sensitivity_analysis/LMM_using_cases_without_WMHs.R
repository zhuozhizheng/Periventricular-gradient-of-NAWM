library(dplyr)
library(lme4)
library(multcomp)
library(emmeans)
library(gplots)
library(lmerTest)
library(ggpubr)

rm(list=ls())
datapath='../data/NODDI_ring_values';
savepath='../Results/NDI';
#savepath='../Results/ODI';

setwd(datapath);

#read the NAWM diffusion data
MRI<-read.csv('NDI_data.csv',header= TRUE);Str<-'NDI'
#MRI<-read.csv('NDI_data.csv',header= TRUE);Str<-'ODI'

rownames(MRI)<-MRI[,1]

col1=1;col2=10
Var_index<-col1:col2

MRI<-MRI[,Var_index]

#read the clincial data
Clinical<-read.csv('Clinical_data.csv',header = TRUE);
rownames(Clinical)<-Clinical[,1]

data<-MRI
data<-data.frame(data);
Clinical<-data.frame(Clinical);

Select_ind1<-which(Clinical[,2]==1&is.na(Clinical$Lesion_volume)) # select cases without WMH in HC groups

data<-MRI[Select_ind1,] 
Clinical<-Clinical[Select_ind1,]

#processing NAN 
NAN_juge<-data[!complete.cases(data)]

for(i in ncol(data))
{ 
   data[,i]<-impute(data[,i],mean) #using the mean value to process the missing values
 
}


#delete the extreme value of each ring in each group
HC_index<-NULL

for (i in 1:ncol(data))
{tem_HC<-data[Clinical[,1]=='1',i]
 

 out<-boxplot.stats(tem_HC)$out
 out_id<-which(tem_HC %in% c(out))
 HC_index<-c(HC_index,out_id)

}

tem_HC<-data[Clinical[,1]=='1',]


index_del<-c(rownames(tem_HC)[HC_index])

index_del<-unique(index_del)

if(length(index_del)>0)
{
data<-data[-which(rownames(data) %in% index_del),]
Clinical<-Clinical[-which(rownames(Clinical) %in% index_del),]
}

#change the path to savepath
setwd(savepath);

#extract clinical and whole brian level data
Group<-data.frame(Clinical[,2])
Age<-Clinical$Age
Sex<-Clinical$Sex..F.1.M.0.
Protocol<-Clinical$Protocol.

MMSE<-Clinical$MMSE
MOCA<-Clinical$MOCA

GM<-Clinical$GM
WM<-Clinical$WM
TIV<-Clinical$TIV

WM_ICVF<-Clinical$WM_ICVF
WM_OD<-Clinical$WM_ODI


Lesion_vol<-Clinical$Lesion_volume


#new_group label
New_group<-data.frame(Group)
New_group[Age<45&Group=='1',]<-'Youth'
New_group[Age>=45&Age<60&Group=='1',]<-'Middle'
New_group[Age>=60&Group=='1',]<-'Older'

data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older'))

#raw parameter plot
ring_num<-length(Var_index)
#prepare data for LMM
New_data0<-data.frame()
for (i in 1:nrow(data))
{for (j in 1:ring_num)
  
{
  New_data0[(i-1)*ring_num+j,1]<-data[i,j]
  New_data0[(i-1)*ring_num+j,2]<-j*3   # 3mm for each ring
  New_data0[(i-1)*ring_num+j,3]<-i
  New_data0[(i-1)*ring_num+j,4]<-New_group[i,1]
  New_data0[(i-1)*ring_num+j,5]<-Age[i]
  New_data0[(i-1)*ring_num+j,6]<-Sex[i]
  New_data0[(i-1)*ring_num+j,7]<-Protocol[i]
  
 
}
  
}

New_data0[,1]<-as.numeric(New_data0[,1])
New_data0[,2]<-as.numeric(New_data0[,2])
New_data0[,3]<-as.numeric(New_data0[,3])
New_data0[,4]<-as.factor(New_data0[,4])
New_data0[,5]<-as.numeric(New_data0[,5])
New_data0[,6]<-as.factor(New_data0[,6])
New_data0[,7]<-as.factor(New_data0[,7])

colnames(New_data0)<-c('value','ring_dist','individual','Group','Age','Sex','Protocol')

pdf(paste0(Str,'_raw_data.pdf'),width=14, height=3)
p<-ggplot(New_data0, aes(Group, value, group=factor(ring_dist), fill=Group)) +
  geom_boxplot() + 
  facet_grid(.~Group, scales = "free_x")+
  scale_y_continuous(limits = c(0.25, 0.8))+theme_bw()+
  theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
print(p)
dev.off()





#z-score transform
data_zscore<-data
for (i in 1:ncol(data))
{
  tem_col<-data[New_group=='Youth',i]

  data_zscore[,i]<-(as.numeric(unlist(data[,i]))-mean(tem_col))/sd(tem_col)

}


Group<-data.frame(Clinical[,2])
Age<-Clinical$Age
Sex<-Clinical$Sex..F.1.M.0.
Protocol<-Clinical$Protocol.

MMSE<-Clinical$MMSE
MOCA<-Clinical$MOCA

GM<-Clinical$GM
WM<-Clinical$WM
TIV<-Clinical$TIV

WM_ICVF<-Clinical$WM_ICVF
WM_OD<-Clinical$WM_ODI

Disease_Duration_years<-Clinical$Disease_Duration_years
Education_years<-Clinical$Education_years

Clinical[,'WM_NODDI_whole_lelvel']<-Clinical$WM_ICVF  
#Clinical[,'WM_NODDI_whole_lelvel']<-Clinical$WM_ODI

Lesion_vol<-Clinical$Lesion_volume


New_group<-data.frame(Group)
New_group[Age<45&Group=='1',]<-'Youth'
New_group[Age>=45&Age<60&Group=='1',]<-'Middle'
New_group[Age>=60&Group=='1',]<-'Elderly'


ring_num<-length(Var_index)
#prepare data for LMM
New_data<-data.frame()
for (i in 1:nrow(data_zscore))
{for (j in 1:ring_num)
  
{
  New_data[(i-1)*ring_num+j,1]<-data_zscore[i,j]
  New_data[(i-1)*ring_num+j,2]<-j*3   # 3mm for each ring
  New_data[(i-1)*ring_num+j,3]<-i
  New_data[(i-1)*ring_num+j,4]<-New_group[i,1]
  New_data[(i-1)*ring_num+j,5]<-Age[i]
  New_data[(i-1)*ring_num+j,6]<-Sex[i]
  New_data[(i-1)*ring_num+j,7]<-Protocol[i]
}
  
}

New_data[,1]<-as.numeric(New_data[,1])
New_data[,2]<-as.numeric(New_data[,2])
New_data[,3]<-as.numeric(New_data[,3])
New_data[,4]<-as.factor(New_data[,4])
New_data[,5]<-as.numeric(New_data[,5])
New_data[,6]<-as.factor(New_data[,6])
New_data[,7]<-as.factor(New_data[,7])


colnames(New_data)<-c('value','ring_dist','individual','Group','Age','Sex','Protocol')
#write.csv(New_data,file = paste0(Str,"_Newdata_results.csv"),
#         quote = F,row.names = F)


#plot before the covariate regressed
New_data$Group<-factor(New_data$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))

pdf(paste0(Str,'_zscore_before_cov_regress.pdf'),width=14, height=3)
p<-ggplot(New_data, aes(Group, value, group=factor(ring_dist), fill=Group)) +
  geom_boxplot() + 
  facet_grid(.~Group, scales = "free_x")+
  scale_y_continuous(limits = c(-3, 3))+theme_bw()+
  theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))

print(p);
dev.off()


LM1<-lm(value~Age+Sex+Protocol,data=New_data)
summary(LM1)

New_data1<-New_data
New_data1$value<-LM1$residuals+mean(New_data1$value)
#plot after the covariate regressed
New_data1$Group<-factor(New_data1$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(paste0(Str,'_zscore_after_cov_regress.pdf'),width=14, height=3)
p<-ggplot(New_data1, aes(Group, value, group=factor(ring_dist), fill=Group)) +
  geom_boxplot() + 
  facet_grid(.~Group, scales = "free_x")+
  scale_y_continuous(limits = c(-3, 3))+theme_bw()+
  theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
print(p);
dev.off()



LME<-lmer(value~Age+Sex+Protocol+ring_dist+Group+ring_dist*Group+(1+ring_dist|individual)+(1+ring_dist|Protocol),data=New_data)


#compared to Youth
Result1<-summary(LME)
#Compared to other Groups
Result2<-summary(glht(LME, linfct=c(
                           'ring_dist:GroupOlder-ring_dist:GroupMiddle=0',
                           'ring_dist:GroupAD-ring_dist:GroupMiddle=0',
                           'ring_dist:GroupPD-ring_dist:GroupMiddle=0',
                           'ring_dist:GroupSVD-ring_dist:GroupMiddle=0',
                           'ring_dist:GroupMS-ring_dist:GroupMiddle=0'
                           
                           )))


Result3<-fix_coef<-fixed.effects(LME)
Result4<-random.effects(LME)

#write.csv(Result3,paste0(Str,'_Result3.csv'))
#write.csv(Result4,paste0(Str,'Result4.csv'))


#prepare the data for forestplot
fit.result<-summary(LME)
df1<-fit.result$coefficients
df2<-confint(LME)
df3<-cbind(df1,df2[5:nrow(df2),])
df4<-data.frame(df3[,c(1,5,6,7)])
df4$Var<-rownames(df4)
colnames(df4)<-c("Coef","Pvalue","Coef_1","Coef_2",'Var')
df5<-df4[,c(5,1,2,3,4)]
df5$Coef_mean<-df5$Coef
df5$Coef<-paste0(round(df5$Coef,2),
               "(",
               round(df5$Coef_1,2),
               "~",
               round(df5$Coef_2,2),
               ")")
df5$Pvalue<-round(df5$Pvalue,3)


write.csv(df5,file =paste0(Str, "_forestplot_results.csv"),
          quote = F,row.names = F)


