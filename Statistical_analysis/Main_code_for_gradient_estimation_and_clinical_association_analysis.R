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

#processing NAN 
NAN_juge<-data[!complete.cases(data)]

for(i in 1:ncol(data))
{ 
   data[,i]<-impute(data[,i],mean) #using the mean value to process the missing values
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
New_group[Group=='3',]<-'AD'
New_group[Group=='4',]<-'PD'
New_group[Group=='5',]<-'SVD'
New_group[Group=='6',]<-'MS'

#clinical variable results summary and plot

#sex
Sex<-as.factor(Sex)
data_cli<-cbind(Sex,New_group)
colnames(data_cli)<-c('Sex','New_group')
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))

pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggplot(data_cli,mapping=aes(x =New_group,fill=Sex,palette = "jco",legend="none"))+
     geom_bar(position="dodge",width=0.75)
print(p1)
dev.off()


My_comparison<-list(c('Middle','Youth'),
                    c('Older','Youth'),
                    c('AD','Youth'),
                    c('PD','Youth'),
                    c('SVD','Youth'),
                    c('MS','Youth'),
                                        
                    c('Older','Middle'),
                    c('AD','Middle'),
                    c('PD','Middle'),
                    c('SVD','Middle'),
                    c('MS','Middle'),
                                        
                    c('AD','Older'),
                    c('PD','Older'),
                    c('SVD','Older'),
                    c('MS','Older'),
                                        
                    c('PD','AD'),
                    c('SVD','AD'),
                    c('MS','AD'),
                                        
                    c('SVD','PD'),
                    c('MS','PD'),
                                        
                    c('MS','SVD'))


#Age
data_cli<-cbind(Age,New_group)

colnames(data_cli)<-c('Age','New_group')

data_cli[,'Age']<-as.numeric(data_cli[,'Age'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))

pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()


#MMSE
data_cli<-cbind(MMSE,New_group)
colnames(data_cli)<-c('MMSE','New_group')
data_cli[,'MMSE']<-as.numeric(data_cli[,'MMSE'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()


#MOCA
data_cli<-cbind(MOCA,New_group)
colnames(data_cli)<-c('MOCA','New_group')
data_cli[,'MOCA']<-as.numeric(data_cli[,'MOCA'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()

#GMV
data_cli<-cbind(GM,New_group)
colnames(data_cli)<-c('GM','New_group')
data_cli[,'GM']<-as.numeric(data_cli[,'GM'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons =My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()

#WMV
data_cli<-cbind(WM,New_group)
colnames(data_cli)<-c('WM','New_group')
data_cli[,'WM']<-as.numeric(data_cli[,'WM'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()

#TIV
data_cli<-cbind(TIV,New_group)
colnames(data_cli)<-c('TIV','New_group')
data_cli[,'TIV']<-as.numeric(data_cli[,'TIV'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()

#WM_ICVF
data_cli<-cbind(WM_ICVF,New_group)
colnames(data_cli)<-c('WM_ICVF','New_group')
data_cli[,'WM_ICVF']<-as.numeric(data_cli[,'WM_ICVF'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()

#WM_OD
data_cli<-cbind(WM_OD,New_group)
colnames(data_cli)<-c('WM_OD','New_group')
data_cli[,'WM_OD']<-as.numeric(data_cli[,'WM_OD'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()

#WMH_vol
data_cli<-cbind(Lesion_vol,New_group)
colnames(data_cli)<-c('Lesion_vol','New_group')
data_cli[,'Lesion_vol']<-as.numeric(data_cli[,'Lesion_vol'])
data_cli$New_group<-factor(data_cli$New_group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(file=paste0(Str,colnames(data_cli)[1],'.pdf'))
p1 = ggboxplot(data_cli, x = "New_group", y = colnames(data_cli)[1],
               color = "New_group", palette = "jco",add = "jitter",legend="none")+
  
  stat_compare_means(comparisons = My_comparison,
                     label="P",
                     p.adjust.method = "BHR",
                     hide.ns=TRUE)
#label.y=c())
#stat_compare_means(method='anova')

print(p1)
dev.off()


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


#delete the extreme value of each ring in each group
HC_index<-NULL
AD_index<-NULL
PD_index<-NULL
SVD_index<-NULL
MS_index<-NULL

#data<-as.numeric(data)

for (i in 1:ncol(data))
{tem_HC<-data_zscore[Clinical[,1]=='1',i]
tem_AD<-data_zscore[Clinical[,1]=='3',i]
tem_PD<-data_zscore[Clinical[,1]=='4',i]
tem_SVD<-data_zscore[Clinical[,1]=='5',i]
tem_MS<-data_zscore[Clinical[,1]=='6',i]

out<-boxplot.stats(tem_HC)$out
out_id<-which(tem_HC %in% c(out))
HC_index<-c(HC_index,out_id)

out<-boxplot.stats(tem_AD)$out
out_id<-which(tem_AD %in% c(out))
AD_index<-c(AD_index,out_id)

out<-boxplot.stats(tem_PD)$out
out_id<-which(tem_PD %in% c(out))
PD_index<-c(PD_index,out_id)

out<-boxplot.stats(tem_SVD)$out
out_id<-which(tem_SVD %in% c(out))
SVD_index<-c(SVD_index,out_id)

out<-boxplot.stats(tem_MS)$out
out_id<-which(tem_MS %in% c(out))
MS_index<-c(MS_index,out_id)

}

tem_HC<-data_zscore[Clinical[,1]=='1',]
tem_AD<-data_zscore[Clinical[,1]=='3',]
tem_PD<-data_zscore[Clinical[,1]=='4',]
tem_SVD<-data_zscore[Clinical[,1]=='5',]
tem_MS<-data_zscore[Clinical[,1]=='6',]

index_del<-c(rownames(tem_HC)[HC_index],rownames(tem_AD)[AD_index],rownames(tem_PD)[PD_index],
             rownames(tem_SVD)[SVD_index],rownames(tem_MS)[MS_index])

index_del<-unique(index_del)

if(length(index_del)>0)
{
  data_zscore<-data_zscore[-which(rownames(data_zscore) %in% index_del),]
  Clinical<-Clinical[-which(rownames(Clinical) %in% index_del),]
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
New_group[Group=='3',]<-'AD'
New_group[Group=='4',]<-'PD'
New_group[Group=='5',]<-'SVD'
New_group[Group=='6',]<-'MS'


#prepare data for LMM
ring_num<-length(Var_index)

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
                           'ring_dist:GroupMS-ring_dist:GroupMiddle=0',
                           
                           
                           'ring_dist:GroupAD-ring_dist:GroupOlder=0',
                           'ring_dist:GroupPD-ring_dist:GroupOlder=0',
                           'ring_dist:GroupSVD-ring_dist:GroupOlder=0',
                           'ring_dist:GroupMS-ring_dist:GroupOlder=0',
                           
                           'ring_dist:GroupPD-ring_dist:GroupAD=0',
                           'ring_dist:GroupSVD-ring_dist:GroupAD=0',
                           'ring_dist:GroupMS-ring_dist:GroupAD=0',
                           
                           'ring_dist:GroupSVD-ring_dist:GroupPD=0',
                           'ring_dist:GroupMS-ring_dist:GroupPD=0',
                           
                           'ring_dist:GroupMS-ring_dist:GroupSVD=0'
                          
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




#clincial correlation analysis in maxin txt
Cor_Group<-New_group
Cor_data<-cbind(Result4$individual$`(Intercept)`,Result4$individual$ring_dist,Age,Sex,TIV,GM,WM,Lesion_vol,
                MMSE,MOCA,Cor_Group,Protocol,Disease_Duration_years,Education_years,WM_NODDI_whole_lelvel)
colnames(Cor_data)<-c('Intercept','Slope','Age','Sex','TIV','GMV','WMV','Lesion_volume',
                      'MMSE','MOCA','Group','Protocol','Disease_Duration','Education','WM_NODDI_whole_lelvel')

Cor_data$Sex<-as.factor(Cor_data$Sex)
Cor_data$Group<-as.factor(Cor_data$Group)
Cor_data$Protocol<-as.factor(Cor_data$Protocol)


#Age with slope

Age_results<-data.frame()

fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Youth'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[1,1]<-R1$coefficients[2,1]
Age_results[1,2]<-R2[2,1]
Age_results[1,3]<-R2[2,2]
Age_results[1,4]<-R1$coefficients[2,3]
Age_results[1,5]<-R1$coefficients[2,4]
Age_results[1,6]<-R1$adj.r.squared

fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Middle'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[2,1]<-R1$coefficients[2,1]
Age_results[2,2]<-R2[2,1]
Age_results[2,3]<-R2[2,2]
Age_results[2,4]<-R1$coefficients[2,3]
Age_results[2,5]<-R1$coefficients[2,4]
Age_results[2,6]<-R1$adj.r.squared


fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Older'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[3,1]<-R1$coefficients[2,1]
Age_results[3,2]<-R2[2,1]
Age_results[3,3]<-R2[2,2]
Age_results[3,4]<-R1$coefficients[2,3]
Age_results[3,5]<-R1$coefficients[2,4]
Age_results[3,6]<-R1$adj.r.squared


fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='AD'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[4,1]<-R1$coefficients[2,1]
Age_results[4,2]<-R2[2,1]
Age_results[4,3]<-R2[2,2]
Age_results[4,4]<-R1$coefficients[2,3]
Age_results[4,5]<-R1$coefficients[2,4]
Age_results[4,6]<-R1$adj.r.squared


fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='PD'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[5,1]<-R1$coefficients[2,1]
Age_results[5,2]<-R2[2,1]
Age_results[5,3]<-R2[2,2]
Age_results[5,4]<-R1$coefficients[2,3]
Age_results[5,5]<-R1$coefficients[2,4]
Age_results[5,6]<-R1$adj.r.squared


fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='SVD'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[6,1]<-R1$coefficients[2,1]
Age_results[6,2]<-R2[2,1]
Age_results[6,3]<-R2[2,2]
Age_results[6,4]<-R1$coefficients[2,3]
Age_results[6,5]<-R1$coefficients[2,4]
Age_results[6,6]<-R1$adj.r.squared


fit0<-lm(Slope~Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='MS'))
R1<-summary(fit0)
R2<-confint(fit0)

Age_results[7,1]<-R1$coefficients[2,1]
Age_results[7,2]<-R2[2,1]
Age_results[7,3]<-R2[2,2]
Age_results[7,4]<-R1$coefficients[2,3]
Age_results[7,5]<-R1$coefficients[2,4]
Age_results[7,6]<-R1$adj.r.squared


write.csv(Age_results,file =paste0(Str, "_Age_results.csv"),
          quote = F,row.names = F)

Cor_data$Group<-factor(Cor_data$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))

pdf(paste0(Str,'_age_association.pdf'),width=14, height=3)
p <- ggplot(Cor_data, aes(Age, Slope)) +
  geom_point(aes(color = Group), size = 1, alpha = 1) +
  facet_wrap(~Group,scales = "free_x",nrow = 1,ncol = 7)+theme_bw()+
  scale_y_continuous(limits = c(-0.05, 0.1))+
  stat_smooth( aes(color = Group, fill = Group), method = "gam",formula =y~x )
print(p)
dev.off()


#GMV with slope
GMV_results<-data.frame()

fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Youth'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[1,1]<-R1$coefficients[2,1]
GMV_results[1,2]<-R2[2,1]
GMV_results[1,3]<-R2[2,2]
GMV_results[1,4]<-R1$coefficients[2,3]
GMV_results[1,5]<-R1$coefficients[2,4]
GMV_results[1,6]<-R1$adj.r.squared

fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Middle'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[2,1]<-R1$coefficients[2,1]
GMV_results[2,2]<-R2[2,1]
GMV_results[2,3]<-R2[2,2]
GMV_results[2,4]<-R1$coefficients[2,3]
GMV_results[2,5]<-R1$coefficients[2,4]
GMV_results[2,6]<-R1$adj.r.squared


fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Older'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[3,1]<-R1$coefficients[2,1]
GMV_results[3,2]<-R2[2,1]
GMV_results[3,3]<-R2[2,2]
GMV_results[3,4]<-R1$coefficients[2,3]
GMV_results[3,5]<-R1$coefficients[2,4]
GMV_results[3,6]<-R1$adj.r.squared


fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='AD'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[4,1]<-R1$coefficients[2,1]
GMV_results[4,2]<-R2[2,1]
GMV_results[4,3]<-R2[2,2]
GMV_results[4,4]<-R1$coefficients[2,3]
GMV_results[4,5]<-R1$coefficients[2,4]
GMV_results[4,6]<-R1$adj.r.squared


fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='PD'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[5,1]<-R1$coefficients[2,1]
GMV_results[5,2]<-R2[2,1]
GMV_results[5,3]<-R2[2,2]
GMV_results[5,4]<-R1$coefficients[2,3]
GMV_results[5,5]<-R1$coefficients[2,4]
GMV_results[5,6]<-R1$adj.r.squared


fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='SVD'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[6,1]<-R1$coefficients[2,1]
GMV_results[6,2]<-R2[2,1]
GMV_results[6,3]<-R2[2,2]
GMV_results[6,4]<-R1$coefficients[2,3]
GMV_results[6,5]<-R1$coefficients[2,4]
GMV_results[6,6]<-R1$adj.r.squared


fit0<-lm(Slope~GMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='MS'))
R1<-summary(fit0)
R2<-confint(fit0)

GMV_results[7,1]<-R1$coefficients[2,1]
GMV_results[7,2]<-R2[2,1]
GMV_results[7,3]<-R2[2,2]
GMV_results[7,4]<-R1$coefficients[2,3]
GMV_results[7,5]<-R1$coefficients[2,4]
GMV_results[7,6]<-R1$adj.r.squared


write.csv(GMV_results,file =paste0(Str, "_GMV_results.csv"),
          quote = F,row.names = F)
Cor_data$Group<-factor(Cor_data$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(paste0(Str,'_GMV_association.pdf'),width=14, height=3)
p <- ggplot(Cor_data, aes(GMV, Slope)) +
  geom_point(aes(color = Group), size = 1, alpha = 1) +
  facet_wrap(~Group,scales = "free_x",nrow = 1,ncol = 7)+theme_bw()+
  scale_y_continuous(limits = c(-0.05, 0.1))+
  stat_smooth( aes(color = Group, fill = Group), method = "gam",formula =y~x )
print(p)
dev.off()


#WMV with slope
WMV_results<-data.frame()

fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Youth'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[1,1]<-R1$coefficients[2,1]
WMV_results[1,2]<-R2[2,1]
WMV_results[1,3]<-R2[2,2]
WMV_results[1,4]<-R1$coefficients[2,3]
WMV_results[1,5]<-R1$coefficients[2,4]
WMV_results[1,6]<-R1$adj.r.squared

fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Middle'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[2,1]<-R1$coefficients[2,1]
WMV_results[2,2]<-R2[2,1]
WMV_results[2,3]<-R2[2,2]
WMV_results[2,4]<-R1$coefficients[2,3]
WMV_results[2,5]<-R1$coefficients[2,4]
WMV_results[2,6]<-R1$adj.r.squared


fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Older'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[3,1]<-R1$coefficients[2,1]
WMV_results[3,2]<-R2[2,1]
WMV_results[3,3]<-R2[2,2]
WMV_results[3,4]<-R1$coefficients[2,3]
WMV_results[3,5]<-R1$coefficients[2,4]
WMV_results[3,6]<-R1$adj.r.squared


fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='AD'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[4,1]<-R1$coefficients[2,1]
WMV_results[4,2]<-R2[2,1]
WMV_results[4,3]<-R2[2,2]
WMV_results[4,4]<-R1$coefficients[2,3]
WMV_results[4,5]<-R1$coefficients[2,4]
WMV_results[4,6]<-R1$adj.r.squared


fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='PD'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[5,1]<-R1$coefficients[2,1]
WMV_results[5,2]<-R2[2,1]
WMV_results[5,3]<-R2[2,2]
WMV_results[5,4]<-R1$coefficients[2,3]
WMV_results[5,5]<-R1$coefficients[2,4]
WMV_results[5,6]<-R1$adj.r.squared


fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='SVD'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[6,1]<-R1$coefficients[2,1]
WMV_results[6,2]<-R2[2,1]
WMV_results[6,3]<-R2[2,2]
WMV_results[6,4]<-R1$coefficients[2,3]
WMV_results[6,5]<-R1$coefficients[2,4]
WMV_results[6,6]<-R1$adj.r.squared


fit0<-lm(Slope~WMV+Age+Sex+TIV+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='MS'))
R1<-summary(fit0)
R2<-confint(fit0)

WMV_results[7,1]<-R1$coefficients[2,1]
WMV_results[7,2]<-R2[2,1]
WMV_results[7,3]<-R2[2,2]
WMV_results[7,4]<-R1$coefficients[2,3]
WMV_results[7,5]<-R1$coefficients[2,4]
WMV_results[7,6]<-R1$adj.r.squared



write.csv(WMV_results,file =paste0(Str, "_WMV_results.csv"),
          quote = F,row.names = F)
Cor_data$Group<-factor(Cor_data$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(paste0(Str,'_WMV_association.pdf'),width=14, height=3)
p <- ggplot(Cor_data, aes(WMV, Slope)) +
  geom_point(aes(color = Group), size = 1, alpha = 1) +
  facet_wrap(~Group,scales = "free_x",nrow = 1,ncol = 7)+theme_bw()+
  scale_y_continuous(limits = c(-0.05, 0.1))+
  stat_smooth( aes(color = Group, fill = Group), method = "gam",formula =y~x )
print(p)
dev.off()


#Lesion_volume with slope
Lesion_volume_results<-data.frame()
Cor_data$Lesion_volume<-log(Cor_data$Lesion_volume)

fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Youth'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[1,1]<-R1$coefficients[2,1]
Lesion_volume_results[1,2]<-R2[2,1]
Lesion_volume_results[1,3]<-R2[2,2]
Lesion_volume_results[1,4]<-R1$coefficients[2,3]
Lesion_volume_results[1,5]<-R1$coefficients[2,4]
Lesion_volume_results[1,6]<-R1$adj.r.squared

fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Middle'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[2,1]<-R1$coefficients[2,1]
Lesion_volume_results[2,2]<-R2[2,1]
Lesion_volume_results[2,3]<-R2[2,2]
Lesion_volume_results[2,4]<-R1$coefficients[2,3]
Lesion_volume_results[2,5]<-R1$coefficients[2,4]
Lesion_volume_results[2,6]<-R1$adj.r.squared


fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Older'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[3,1]<-R1$coefficients[2,1]
Lesion_volume_results[3,2]<-R2[2,1]
Lesion_volume_results[3,3]<-R2[2,2]
Lesion_volume_results[3,4]<-R1$coefficients[2,3]
Lesion_volume_results[3,5]<-R1$coefficients[2,4]
Lesion_volume_results[3,6]<-R1$adj.r.squared


fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='AD'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[4,1]<-R1$coefficients[2,1]
Lesion_volume_results[4,2]<-R2[2,1]
Lesion_volume_results[4,3]<-R2[2,2]
Lesion_volume_results[4,4]<-R1$coefficients[2,3]
Lesion_volume_results[4,5]<-R1$coefficients[2,4]
Lesion_volume_results[4,6]<-R1$adj.r.squared


fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='PD'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[5,1]<-R1$coefficients[2,1]
Lesion_volume_results[5,2]<-R2[2,1]
Lesion_volume_results[5,3]<-R2[2,2]
Lesion_volume_results[5,4]<-R1$coefficients[2,3]
Lesion_volume_results[5,5]<-R1$coefficients[2,4]
Lesion_volume_results[5,6]<-R1$adj.r.squared


fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='SVD'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[6,1]<-R1$coefficients[2,1]
Lesion_volume_results[6,2]<-R2[2,1]
Lesion_volume_results[6,3]<-R2[2,2]
Lesion_volume_results[6,4]<-R1$coefficients[2,3]
Lesion_volume_results[6,5]<-R1$coefficients[2,4]
Lesion_volume_results[6,6]<-R1$adj.r.squared


fit0<-lm(Slope~Lesion_volume+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='MS'))
R1<-summary(fit0)
R2<-confint(fit0)

Lesion_volume_results[7,1]<-R1$coefficients[2,1]
Lesion_volume_results[7,2]<-R2[2,1]
Lesion_volume_results[7,3]<-R2[2,2]
Lesion_volume_results[7,4]<-R1$coefficients[2,3]
Lesion_volume_results[7,5]<-R1$coefficients[2,4]
Lesion_volume_results[7,6]<-R1$adj.r.squared


write.csv(Lesion_volume_results,file =paste0(Str, "_Lesion_volume_results.csv"),
          quote = F,row.names = F)
Cor_data$Group<-factor(Cor_data$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(paste0(Str,'_Lesion_volume_association.pdf'),width=14, height=3)
p <- ggplot(Cor_data, aes(Lesion_volume, Slope)) +
  geom_point(aes(color = Group), size = 1, alpha = 1) +
  facet_wrap(~Group,scales = "free_x",nrow = 1,ncol = 7)+theme_bw()+
  scale_y_continuous(limits = c(-0.05, 0.1))+
  stat_smooth( aes(color = Group, fill = Group), method = "gam",formula =y~x )
print(p)
dev.off()


#MOCA with slope

MOCA_results<-data.frame()

fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Youth'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[1,1]<-R1$coefficients[2,1]
MOCA_results[1,2]<-R2[2,1]
MOCA_results[1,3]<-R2[2,2]
MOCA_results[1,4]<-R1$coefficients[2,3]
MOCA_results[1,5]<-R1$coefficients[2,4]
MOCA_results[1,6]<-R1$adj.r.squared

fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Middle'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[2,1]<-R1$coefficients[2,1]
MOCA_results[2,2]<-R2[2,1]
MOCA_results[2,3]<-R2[2,2]
MOCA_results[2,4]<-R1$coefficients[2,3]
MOCA_results[2,5]<-R1$coefficients[2,4]
MOCA_results[2,6]<-R1$adj.r.squared


fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Older'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[3,1]<-R1$coefficients[2,1]
MOCA_results[3,2]<-R2[2,1]
MOCA_results[3,3]<-R2[2,2]
MOCA_results[3,4]<-R1$coefficients[2,3]
MOCA_results[3,5]<-R1$coefficients[2,4]
MOCA_results[3,6]<-R1$adj.r.squared


fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='AD'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[4,1]<-R1$coefficients[2,1]
MOCA_results[4,2]<-R2[2,1]
MOCA_results[4,3]<-R2[2,2]
MOCA_results[4,4]<-R1$coefficients[2,3]
MOCA_results[4,5]<-R1$coefficients[2,4]
MOCA_results[4,6]<-R1$adj.r.squared


fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='PD'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[5,1]<-R1$coefficients[2,1]
MOCA_results[5,2]<-R2[2,1]
MOCA_results[5,3]<-R2[2,2]
MOCA_results[5,4]<-R1$coefficients[2,3]
MOCA_results[5,5]<-R1$coefficients[2,4]
MOCA_results[5,6]<-R1$adj.r.squared


fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='SVD'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[6,1]<-R1$coefficients[2,1]
MOCA_results[6,2]<-R2[2,1]
MOCA_results[6,3]<-R2[2,2]
MOCA_results[6,4]<-R1$coefficients[2,3]
MOCA_results[6,5]<-R1$coefficients[2,4]
MOCA_results[6,6]<-R1$adj.r.squared


fit0<-lm(Slope~MOCA+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='MS'))
R1<-summary(fit0)
R2<-confint(fit0)

MOCA_results[7,1]<-R1$coefficients[2,1]
MOCA_results[7,2]<-R2[2,1]
MOCA_results[7,3]<-R2[2,2]
MOCA_results[7,4]<-R1$coefficients[2,3]
MOCA_results[7,5]<-R1$coefficients[2,4]
MOCA_results[7,6]<-R1$adj.r.squared


write.csv(MOCA_results,file =paste0(Str, "_MOCA_results.csv"),
          quote = F,row.names = F)
Cor_data$Group<-factor(Cor_data$Group,levels=c('Youth','Middle','older','AD','PD','SVD','MS'))
pdf(paste0(Str,'_MOCA_association.pdf'),width=14, height=3)
p <- ggplot(Cor_data, aes(MOCA, Slope)) +
  geom_point(aes(color = Group), size = 1, alpha = 1) +
  facet_wrap(~Group,scales = "free_x",nrow = 1,ncol = 7)+theme_bw()+
  scale_y_continuous(limits = c(-0.05, 0.1))+
  stat_smooth( aes(color = Group, fill = Group), method = "gam",formula =y~x )
print(p)
dev.off()



#MMSE with slope

MMSE_results<-data.frame()

fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Youth'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[1,1]<-R1$coefficients[2,1]
MMSE_results[1,2]<-R2[2,1]
MMSE_results[1,3]<-R2[2,2]
MMSE_results[1,4]<-R1$coefficients[2,3]
MMSE_results[1,5]<-R1$coefficients[2,4]
MMSE_results[1,6]<-R1$adj.r.squared

fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Middle'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[2,1]<-R1$coefficients[2,1]
MMSE_results[2,2]<-R2[2,1]
MMSE_results[2,3]<-R2[2,2]
MMSE_results[2,4]<-R1$coefficients[2,3]
MMSE_results[2,5]<-R1$coefficients[2,4]
MMSE_results[2,6]<-R1$adj.r.squared


fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='Older'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[3,1]<-R1$coefficients[2,1]
MMSE_results[3,2]<-R2[2,1]
MMSE_results[3,3]<-R2[2,2]
MMSE_results[3,4]<-R1$coefficients[2,3]
MMSE_results[3,5]<-R1$coefficients[2,4]
MMSE_results[3,6]<-R1$adj.r.squared


fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='AD'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[4,1]<-R1$coefficients[2,1]
MMSE_results[4,2]<-R2[2,1]
MMSE_results[4,3]<-R2[2,2]
MMSE_results[4,4]<-R1$coefficients[2,3]
MMSE_results[4,5]<-R1$coefficients[2,4]
MMSE_results[4,6]<-R1$adj.r.squared


fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='PD'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[5,1]<-R1$coefficients[2,1]
MMSE_results[5,2]<-R2[2,1]
MMSE_results[5,3]<-R2[2,2]
MMSE_results[5,4]<-R1$coefficients[2,3]
MMSE_results[5,5]<-R1$coefficients[2,4]
MMSE_results[5,6]<-R1$adj.r.squared


fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='SVD'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[6,1]<-R1$coefficients[2,1]
MMSE_results[6,2]<-R2[2,1]
MMSE_results[6,3]<-R2[2,2]
MMSE_results[6,4]<-R1$coefficients[2,3]
MMSE_results[6,5]<-R1$coefficients[2,4]
MMSE_results[6,6]<-R1$adj.r.squared


fit0<-lm(Slope~MMSE+Age+Sex+Protocol+WM_NODDI_whole_lelvel,data=Cor_data,subset=(Cor_data$Group=='MS'))
R1<-summary(fit0)
R2<-confint(fit0)

MMSE_results[7,1]<-R1$coefficients[2,1]
MMSE_results[7,2]<-R2[2,1]
MMSE_results[7,3]<-R2[2,2]
MMSE_results[7,4]<-R1$coefficients[2,3]
MMSE_results[7,5]<-R1$coefficients[2,4]
MMSE_results[7,6]<-R1$adj.r.squared



write.csv(MMSE_results,file =paste0(Str, "_MMSE_results.csv"),
          quote = F,row.names = F)
Cor_data$Group<-factor(Cor_data$Group,levels=c('Youth','Middle','Older','AD','PD','SVD','MS'))
pdf(paste0(Str,'_MMSE_association.pdf'),width=14, height=3)
p <- ggplot(Cor_data, aes(MMSE, Slope)) +
  geom_point(aes(color = Group), size = 1, alpha = 1) +
  facet_wrap(~Group,scales = "free_x",nrow = 1,ncol = 7)+theme_bw()+
  scale_y_continuous(limits = c(-0.05, 0.1))+
  stat_smooth( aes(color = Group, fill = Group), method = "gam",formula =y~x )
print(p)
dev.off()



