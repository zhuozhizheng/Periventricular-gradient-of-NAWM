clear all;clc;

[data1 txt1 raw1]=xlsread('...data\Gene_data\Disease_GWAS_genes.xlsx',2);
[data2 txt2 raw2]=xlsread('...data\Gene_data\Disease_GWAS_genes.xlsx',3);
[data3 txt3 raw3]=xlsread('...data\Gene_data\Disease_GWAS_genes.xlsx',4);
[data4 txt4 raw4]=xlsread('...data\Gene_data\Disease_GWAS_genes.xlsx',1);

% NDI
Str='NDI';
savepath='...\Results\NDI\BP'
[data_pls_p txt_pls_p raw_pls_p]=xlsread('...\Results\NDI\BP\NDI_Pos_ALL_GO.xlsx');
[data_pls_n txt_pls_n raw_pls_n]=xlsread('...\Results\NDI\BP\NDI_Neg_ALL_GO.xlsx');

% ODI
% Str='ODI';
% savepath='...\Results\ODI\BP'
% [data_pls_p txt_pls_p raw_pls_p]=xlsread('...\Results\ODI\BP\ODI_Pos_ALL_GO.xlsx');
% [data_pls_n txt_pls_n raw_pls_n]=xlsread('...\Results\ODI\BP\ODI_Neg_ALL_GO.xlsx');


cd(savepath);
AD_gene=unique(raw1(2:end,2));PD_gene=unique(raw2(2:end,2));SVD_gene=unique(raw3(2:end,2));MS_gene=unique(raw4(2:end,2));
tem=intersect(AD_gene,PD_gene);
tem=intersect(tem,SVD_gene);
tem=intersect(tem,MS_gene);

%intersect(bground,genelist_p);

bground=[AD_gene;PD_gene;SVD_gene;MS_gene];
%bground=bground(2:end);

for i=1:size(raw_pls_p,2)
    if i==1;Str1='HC';
    elseif i==2;Str1='AD';
    elseif i==3;Str1='PD'; 
    elseif i==4;Str1='SVD';      
    elseif i==5;Str1='MS';
    end
clear bground_red_pgenelist_p genelist_n;
genelist_p=raw_pls_p(2:end,i); 
num=0;ind=[];
for j=1:length(genelist_p);
    if find(isnan(genelist_p{j}));
       num=num+1;   
        ind(num)=j;
    end
end
genelist_p(ind)=[];

bground_red_p=intersect(bground,genelist_p);

genelist_n=raw_pls_n(2:end,i);

num=0;ind=[];
for j=1:length(genelist_n);
    if find(isnan(genelist_n{j}));
       num=num+1;   
        ind(num)=j;
    end
end
genelist_n(ind)=[];
bground_red_n=intersect(bground,genelist_n);
%AD
Overlap_AD_p=intersect(genelist_p,AD_gene);
tem_p=hygepdf(length(Overlap_AD_p):length(bground_red_p),length(bground),length(AD_gene),length(bground_red_p));
P_overlap_AD_p=sum(tem_p);

Overlap_AD_n=intersect(genelist_n,AD_gene);
tem_p=hygepdf(length(Overlap_AD_n):length(bground_red_n),length(bground),length(AD_gene),length(bground_red_n));
P_overlap_AD_n=sum(tem_p);
%PD
Overlap_PD_p=intersect(genelist_p,PD_gene);
tem_p=hygepdf(length(Overlap_PD_p):length(bground_red_p),length(bground),length(PD_gene),length(bground_red_p));
P_overlap_PD_p=sum(tem_p);

Overlap_PD_n=intersect(genelist_n,PD_gene);
tem_p=hygepdf(length(Overlap_PD_n):length(bground_red_n),length(bground),length(PD_gene),length(bground_red_n));
P_overlap_PD_n=sum(tem_p);

%SVD
Overlap_SVD_p=intersect(genelist_p,SVD_gene);
tem_p=hygepdf(length(Overlap_SVD_p):length(bground_red_p),length(bground),length(SVD_gene),length(bground_red_p));
P_overlap_SVD_p=sum(tem_p);

Overlap_SVD_n=intersect(genelist_n,SVD_gene);
tem_p=hygepdf(length(Overlap_SVD_n):length(bground_red_n),length(bground),length(SVD_gene),length(bground_red_n));
P_overlap_SVD_n=sum(tem_p);

%MS
Overlap_MS_p=intersect(genelist_p,MS_gene);
tem_p=hygepdf(length(Overlap_MS_p):length(bground_red_p),length(bground),length(MS_gene),length(bground_red_p));
P_overlap_MS_p=sum(tem_p);

Overlap_MS_n=intersect(genelist_n,MS_gene);
tem_p=hygepdf(length(Overlap_MS_n):length(bground_red_n),length(bground),length(MS_gene),length(bground_red_n));
P_overlap_MS_n=sum(tem_p);

Neg_FDR_p_gene_overlap{i}=mafdr([P_overlap_AD_n,P_overlap_PD_n,P_overlap_SVD_n,P_overlap_MS_n],'BHFDR','true');
Pos_FDR_p_gene_overlap{i}=mafdr([P_overlap_AD_p,P_overlap_PD_p,P_overlap_SVD_p,P_overlap_MS_p],'BHFDR','true')
%save(['Disease_overlap_' Str '_' Str1 '.mat']);
end


% Neg_FDR_p_gene_overlap=mafdr([P_overlap_AD_n,P_overlap_PD_n,P_overlap_SVD_n,P_overlap_MS_n],'BHFDR','true');
% Pos_FDR_p_gene_overlap=mafdr([P_overlap_AD_p,P_overlap_PD_p,P_overlap_SVD_p,P_overlap_MS_p],'BHFDR','true')




for i=1:size(raw_pls_p,2)
    if i==1;Str1='HC';
    elseif i==2;Str1='AD';
    elseif i==3;Str1='PD'; 
    elseif i==4;Str1='SVD';      
    elseif i==5;Str1='MS';
    end
clear bground_red_pgenelist_p genelist_n;
genelist_p=raw_pls_p(2:end,i); 
num=0;ind=[];
for j=1:length(genelist_p);
    if find(isnan(genelist_p{j}));
       num=num+1;   
        ind(num)=j;
    end
end
genelist_p(ind)=[];

bground_red_p=intersect(bground,genelist_p);

genelist_n=raw_pls_n(2:end,i);

num=0;ind=[];
for j=1:length(genelist_n);
    if find(isnan(genelist_n{j}));
       num=num+1;   
        ind(num)=j;
    end
end
genelist_n(ind)=[];
bground_red_n=intersect(bground,genelist_n);
%AD
Overlap_AD_p=intersect(genelist_p,AD_gene);
tem_p=hygepdf(length(Overlap_AD_p):length(bground_red_p),length(bground),length(AD_gene),length(bground_red_p));
P_overlap_AD_p=sum(tem_p);

Overlap_AD_n=intersect(genelist_n,AD_gene);
tem_p=hygepdf(length(Overlap_AD_n):length(bground_red_n),length(bground),length(AD_gene),length(bground_red_n));
P_overlap_AD_n=sum(tem_p);
%PD
Overlap_PD_p=intersect(genelist_p,PD_gene);
tem_p=hygepdf(length(Overlap_PD_p):length(bground_red_p),length(bground),length(PD_gene),length(bground_red_p));
P_overlap_PD_p=sum(tem_p);

Overlap_PD_n=intersect(genelist_n,PD_gene);
tem_p=hygepdf(length(Overlap_PD_n):length(bground_red_n),length(bground),length(PD_gene),length(bground_red_n));
P_overlap_PD_n=sum(tem_p);

%SVD
Overlap_SVD_p=intersect(genelist_p,SVD_gene);
tem_p=hygepdf(length(Overlap_SVD_p):length(bground_red_p),length(bground),length(SVD_gene),length(bground_red_p));
P_overlap_SVD_p=sum(tem_p);

Overlap_SVD_n=intersect(genelist_n,SVD_gene);
tem_p=hygepdf(length(Overlap_SVD_n):length(bground_red_n),length(bground),length(SVD_gene),length(bground_red_n));
P_overlap_SVD_n=sum(tem_p);

%MS
Overlap_MS_p=intersect(genelist_p,MS_gene);
tem_p=hygepdf(length(Overlap_MS_p):length(bground_red_p),length(bground),length(MS_gene),length(bground_red_p));
P_overlap_MS_p=sum(tem_p);

Overlap_MS_n=intersect(genelist_n,MS_gene);
tem_p=hygepdf(length(Overlap_MS_n):length(bground_red_n),length(bground),length(MS_gene),length(bground_red_n));
P_overlap_MS_n=sum(tem_p);

Cell_Results_Neg{(i-1)*5+1,1}='AD';
Cell_Results_Neg{(i-1)*5+1,2}=length(AD_gene);
Cell_Results_Neg{(i-1)*5+1,3}=length(genelist_n);
Cell_Results_Neg{(i-1)*5+1,4}=length(Overlap_AD_n);
Cell_Results_Neg{(i-1)*5+1,5}=P_overlap_AD_n;
Cell_Results_Neg{(i-1)*5+1,6}=Neg_FDR_p_gene_overlap{i}(1);
Cell_Results_Neg{(i-1)*5+1,7}=Str1;

Cell_Results_Neg{(i-1)*5+2,1}='PD';
Cell_Results_Neg{(i-1)*5+2,2}=length(PD_gene);
Cell_Results_Neg{(i-1)*5+2,3}=length(genelist_n);
Cell_Results_Neg{(i-1)*5+2,4}=length(Overlap_PD_n);
Cell_Results_Neg{(i-1)*5+2,5}=P_overlap_PD_n;
Cell_Results_Neg{(i-1)*5+2,6}=Neg_FDR_p_gene_overlap{i}(2);
Cell_Results_Neg{(i-1)*5+2,7}=Str1;


Cell_Results_Neg{(i-1)*5+3,1}='SVD';
Cell_Results_Neg{(i-1)*5+3,2}=length(SVD_gene);
Cell_Results_Neg{(i-1)*5+3,3}=length(genelist_n);
Cell_Results_Neg{(i-1)*5+3,4}=length(Overlap_SVD_n);
Cell_Results_Neg{(i-1)*5+3,5}=P_overlap_SVD_n;
Cell_Results_Neg{(i-1)*5+3,6}=Neg_FDR_p_gene_overlap{i}(3);
Cell_Results_Neg{(i-1)*5+3,7}=Str1;

Cell_Results_Neg{(i-1)*5+4,1}='MS';
Cell_Results_Neg{(i-1)*5+4,2}=length(MS_gene);
Cell_Results_Neg{(i-1)*5+4,3}=length(genelist_n);
Cell_Results_Neg{(i-1)*5+4,4}=length(Overlap_MS_n);
Cell_Results_Neg{(i-1)*5+4,5}=P_overlap_MS_n;
Cell_Results_Neg{(i-1)*5+4,6}=Neg_FDR_p_gene_overlap{i}(4);
Cell_Results_Neg{(i-1)*5+4,7}=Str1;


Cell_Results_Pos{(i-1)*5+1,1}='AD';
Cell_Results_Pos{(i-1)*5+1,2}=length(AD_gene);
Cell_Results_Pos{(i-1)*5+1,3}=length(genelist_p);
Cell_Results_Pos{(i-1)*5+1,4}=length(Overlap_AD_p);
Cell_Results_Pos{(i-1)*5+1,5}=P_overlap_AD_p;
Cell_Results_Pos{(i-1)*5+1,6}=Pos_FDR_p_gene_overlap{i}(1);
Cell_Results_Pos{(i-1)*5+1,7}=Str1;

Cell_Results_Pos{(i-1)*5+2,1}='PD';
Cell_Results_Pos{(i-1)*5+2,2}=length(PD_gene);
Cell_Results_Pos{(i-1)*5+2,3}=length(genelist_p);
Cell_Results_Pos{(i-1)*5+2,4}=length(Overlap_PD_p);
Cell_Results_Pos{(i-1)*5+2,5}=P_overlap_PD_p;
Cell_Results_Pos{(i-1)*5+2,6}=Pos_FDR_p_gene_overlap{i}(2);
Cell_Results_Pos{(i-1)*5+2,7}=Str1;


Cell_Results_Pos{(i-1)*5+3,1}='SVD';
Cell_Results_Pos{(i-1)*5+3,2}=length(SVD_gene);
Cell_Results_Pos{(i-1)*5+3,3}=length(genelist_p);
Cell_Results_Pos{(i-1)*5+3,4}=length(Overlap_SVD_p);
Cell_Results_Pos{(i-1)*5+3,5}=P_overlap_SVD_p;
Cell_Results_Pos{(i-1)*5+3,6}=Pos_FDR_p_gene_overlap{i}(3);
Cell_Results_Pos{(i-1)*5+3,7}=Str1;

Cell_Results_Pos{(i-1)*5+4,1}='MS';
Cell_Results_Pos{(i-1)*5+4,2}=length(MS_gene);
Cell_Results_Pos{(i-1)*5+4,3}=length(genelist_p);
Cell_Results_Pos{(i-1)*5+4,4}=length(Overlap_MS_p);
Cell_Results_Pos{(i-1)*5+4,5}=P_overlap_MS_p;
Cell_Results_Pos{(i-1)*5+4,6}=Pos_FDR_p_gene_overlap{i}(4);
Cell_Results_Pos{(i-1)*5+4,7}=Str1;
%save(['Disease_overlap_' Str '_' Str1 '.mat']);
end

savepath='...\Results\NDI'  
cd(savepath);
xlswrite('Neg_disease_result.xlsx',Cell_Results_Neg,'Cell_Results_Neg')
xlswrite('Pos_disease_result.xlsx',Cell_Results_Pos,'Cell_Results_Pos')

