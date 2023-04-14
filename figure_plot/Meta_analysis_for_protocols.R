library(metafor);
library(forestplot)

datapath='D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/meta_analysis'

savepath='D:/ZZZ/Manuscripts/Gradient_study/lifespan_gradient/外部修改/Elsvier/Nature Communication/Solution/Reviewer2/meta_analysis/results'

setwd(datapath);

data<-read.csv('data.csv',header=T);

colnames(data)<-c('Group','Metric','Protocol','LCI','UCI','Coef','Pvalue')
data[,'SE']<-(data[,'UCI']-data[,'LCI'])/(2*1.96)

Stats_result<-data.frame(matrix(0,6*2,12));

setwd(savepath);

for (i1 in 1:6)
{
  i2
  if(i1==2){Group1='Older'};
  if(i1==3){Group1='AD'};
  if(i1==4){Group1='PD'};
  if(i1==5){Group1='SVD'};
  if(i1==6){Group1='MS'};
  
  for (i2 in 1:2)
  {
    if(i2==1){Metric1='NDI'};
    if(i2==2){Metric1='ODI'};
    
    tem_data<-data[data$Group==Group1&data$Metric==Metric1,];
    
    metamod<-rma(yi=Coef,data=tem_data,sei=SE,method='DL')
    
    if (metamod$I2<=50&metamod$QEp>=0.1)
    {
      metamod<-rma(yi=Coef,data=tem_data,sei=SE,method='FE')
      
      
      pdf(paste0(Group1,'_',Metric1,'_meta_analysis.pdf'),width=8, height=4)
      forestplot<-forest(metamod,refline=0,mlab="Random-effect Model",
                         slab=tem_data$Protocol,xlab='Coef',showweights = T)
      title(paste0(Group1,' ',Metric1,' Meta Results'))
     
      dev.off()
    }
    else
    {
      metamod<-rma(yi=Coef,data=tem_data,sei=SE,method='DL')
      
      
      pdf(paste0(Group1,'_',Metric1,'_meta_analysis.pdf'),width=8, height=4)
      forestplot<-forest(metamod,refline=0,mlab="Fixed-effect Model",
                         slab=tem_data$Protocol,xlab='Coef',showweights = T)
      title(paste0(Group1,' ',Metric1,' Meta Results'))
      dev.off()
    }
      
    tem<-summary(metamod)
    
    Stats_result[(i1-1)*2+i2,1]<-Group1;
    Stats_result[(i1-1)*2+i2,2]<-Metric1;
    Stats_result[(i1-1)*2+i2,3]<-tem$QE
    Stats_result[(i1-1)*2+i2,4]<-tem$QEp;
    Stats_result[(i1-1)*2+i2,5]<-tem$I2;
    Stats_result[(i1-1)*2+i2,6]<-tem$tau2;
    Stats_result[(i1-1)*2+i2,7]<-tem$H2;
    Stats_result[(i1-1)*2+i2,8]<-tem$method;
    Stats_result[(i1-1)*2+i2,9]<-tem$beta;
    Stats_result[(i1-1)*2+i2,10]<-tem$ci.lb;
    Stats_result[(i1-1)*2+i2,11]<-tem$ci.ub;
    Stats_result[(i1-1)*2+i2,12]<-tem$pval;
    
  }
}
colnames(Stats_result)<-c('Group','Mertic','Q','Qp',
                         'I2','tau2','H2','Method',
                         'Beta_meta','LCI','UCI','Pvalue')

write.csv(Stats_result,paste0('Meta_summary_for_diseases.csv'));


#according to I^2 <50% and Q p>0.1 using FE, otherwise using DL

#fixed model method='FE'
#random-effect model method='DL'

# metamod<-rma(yi=Coef,data=data,sei=SE,method='DL')
# 
# forestplot<-forest(metamod,refline=1,mlab="Random-effect Model",
#                   slab=data$锘縋rotocol,xlab='Coef',showweights = T)
# text(-500,6:1,pos=2,data$Scanner)
# text(c(-1600,-500,300,800),8,pos=c(4,2,4,4),c('Protocol ID','Scanner','Weights','β[95%CI]'),cex=1,font=2)
# 
# 
# funnel(metamod)
# ranktest(metamod)
# regtest(metamod)
# 
# leavelout(metamod,digits = 3)


