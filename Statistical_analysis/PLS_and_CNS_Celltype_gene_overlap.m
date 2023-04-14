clear all;clc;
% %NDI
datapath='...\data\Gene_data';
cd(datapath);
[data1 txt1 raw1]=xlsread('Perventricular_NDI_for_PLS.csv',1);
savepath='...\Gene_Results\NDI';

% %ODI
% datapath='...\data\Gene_data';
% cd(datapath);
% [data1 txt1 raw1]=xlsread('Perventricular_ODI_for_PLS.csv,1);
% savepath='...\Gene_Results\ODI';


disp('data loading ....');

#[data2 txt2 raw2]=xlsread('gene_expraession_in_perventricular_rings.csv',1);
cell_data_ring=importdata('gene_expraession_in_perventricular_rings.csv');
raw_cell_ring=cell_data_ring.textdata(1,2:end);
genename=raw_cell_ring;

col_start=1;col_end=10;

Y=data1(1,:)';X=data2(col_start:col_end,2:end);Str='HC';
%Y=data1(2,:)';X=data2(col_start:col_end,2:end);Str='AD';
%Y=data1(3,:)';X=data2(col_start:col_end,2:end);Str='PD';
%Y=data1(4,:)';X=data2(col_start:col_end,2:end);Str='SVD';
%Y=data1(5,:)';X=data2(col_start:col_end,2:end);Str='MS';

Y(isnan(X(:,1))==1,:)=[];X(isnan(X(:,1))==1,:)=[];

cd(savepath);
disp('Gene screening by PLS ....');
Y(isnan(X(:,1))==1,:)=[];X(isnan(X(:,1))==1,:)=[];

Y0=zscore(Y);X0=zscore(X,0,1);
[XL YL XS YS BETA PCTVAR MSE stats]=plsregress(X0,Y0);

cd(savepath);

% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:
rep=5000; 
dim=1;
[XL YL XS YS BETA PCTVAR MSE stats]=plsregress(X0,Y0);
temp=cumsum(100*PCTVAR(2,1:dim));

Rsquared = temp(dim);
    parfor i=1:rep
        %i
        order=randperm(size(Y0,1));
        Yp=Y0(order,:);

        [XL1,YL1,XS1,YS1,BETA1,PCTVAR1,MSE1,stats1]=plsregress(X0,Yp,dim);

        temp=cumsum(100*PCTVAR1(2,1:dim));
        Rsq(i) = temp(dim);
    end
    dim
    R(dim)=Rsquared;
    p_PCTVAR(dim)=length(find(Rsq>=Rsquared))/rep;

 
rep=100000;    
Y_ind=randi(length(Y),rep,length(Y));

parfor i=1:rep
   %i  
   [XL1 YL1 XS1 YS1 BETA1 PCTVAR1 MSE1 stats1]=plsregress(X0(Y_ind(i,:),:),Y0(Y_ind(i,:)),dim);
   
   Per_W(:,i)=(stats1.W(:,1));
end


Per_W=[Per_W, stats.W(:,1)];
Z_W=stats.W(:,1)./std(Per_W,0,2);%Z-score calculation
Z_W_5=prctile(Z_W,5,2);
Z_W_95=prctile(Z_W,95,2);


parfor i=1:length(Z_W);
[H_W P_W(i,1)]=ztest(Z_W(i,1),0,1);
end

[FDR_W Q_W]=mafdr(P_W,'Lambda',[0.01:0.01:0.95]);  

Z_weight=stats.W(:,1)./std(XS(:,1));

[Z_W_sort Ind]=sort(Z_W,'ascend');

genename1=genename(Ind);
Z_weight=Z_weight(Ind);
X1=X(:,Ind);

ind_P=find(Z_W_sort>5);% 
ind_N=find(Z_W_sort<-5);% 

P_label=Z_W_sort(ind_P);
N_label=Z_W_sort(ind_N);

P_genename=genename1(ind_P);
N_genename=genename1(ind_N);

P_gene=X1(:,ind_P);
N_gene=X1(:,ind_N);

Genename=[P_genename N_genename];
Gene_list=[P_gene N_gene];
Label=[P_label;N_label];

for i=1:length(Genename);
    list_gene{i,1}=Genename{i};
    list_gene{i,2}=Label(i);
    [list_gene{i,3} list_gene{i,4}]=corr(Gene_list(:,i),Y);
end

cd(savepath);xlswrite([Str '_Gene_list_information.csv'],Gene_list',1);
cd(savepath);xlswrite([Str '_Genename_column.csv'],list_gene,1);


%important information saving
[tem_value Ind_max_cor]=max(cell2mat(list_gene(:,3)));

figure('color','w'),plot(Gene_list(:,Ind_max_cor),Y,'*');
Xz=Gene_list(:,Ind_max_cor);Yz=Y;
[b dev stats]=glmfit(Xz,Yz);
Yhat=glmval(b,Xz,'identity');

[RHO_max_cor PVAL_max_cor]=corr(Gene_list(:,Ind_max_cor),Y);


%permutation test
for i=1:rep 

    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [RHO_tem(i) PVAL_tem]=corr(Gene_list(:,Ind_max_cor),Yp);
end

if RHO_max_cor>0
p_RHO_max_cor=length(find(RHO_max_cor<=RHO_tem))/rep;
else
 p_RHO_max_cor=length(find(RHO_max_cor>=RHO_tem))/rep;   
end
%permutation test

hold on;plot(Gene_list(:,Ind_max_cor),Yhat,'-r');title(['R=',num2str(RHO_max_cor),'--P=',num2str(p_RHO_max_cor)]);
xlabel('Gene-expression');ylabel('T-value');
print(gcf,[Str,'_gene_image_assocaition_max_pos_cor.dng'],'-dpng');
savefig(gcf,[Str,'_gene_image_assocoaition_max_pos_cor.dng.fig']);
close all;

%maping max_pos_gene_expression into volume
tem_imagedata=zeros(size(imagedata));
for i=1:length(Gene_list(:,Ind_max_cor))
    tem_imagedata(imagedata==(i+col_start-1))=Gene_list(i,Ind_max_cor);
end
tem.fname=[Str '_max_pos_gene_expression.nii'];
spm_write_vol(tem,tem_imagedata);
%maping max_pos_gene_expression into volume


[tem_value Ind_min_cor]=min(cell2mat(list_gene(:,3)));

figure('color','w'),plot(Gene_list(:,Ind_min_cor),Y,'*');
Xz=Gene_list(:,Ind_min_cor);Yz=Y;
[b dev stats]=glmfit(Xz,Yz);
Yhat=glmval(b,Xz,'identity');

[RHO_min_cor PVAL_min_cor]=corr(Gene_list(:,Ind_min_cor),Y);

%permutation test
for i=1:rep 

    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [RHO_tem(i) PVAL_tem]=corr(Gene_list(:,Ind_min_cor),Yp);
end
if RHO_min_cor>0
p_RHO_min_cor=length(find(RHO_min_cor<=RHO_tem))/rep;
else
 p_RHO_min_cor=length(find(RHO_min_cor>=RHO_tem))/rep;   
end
%permutation test

hold on;plot(Gene_list(:,Ind_min_cor),Yhat,'-r');title(['R=',num2str(RHO_min_cor),'--P=',num2str(p_RHO_min_cor)]);
xlabel('Gene-expression');ylabel('T-value');
print(gcf,[Str,'_gene_image_assocaition_max_neg_cor.dng'],'-dpng');
savefig(gcf,[Str,'_gene_image_assocoaition_max_neg_cor.dng.fig']);
close all;


[RHO PVAL]=corr(XS(:,1),Y);

figure('color','w'),plot(XS(:,1),Y,'*');
Xz=XS(:,1);Yz=Y;
[b dev stats]=glmfit(Xz,Yz);
Yhat=glmval(b,Xz,'identity');

%permutation test
for i=1:rep 
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [RHO_tem(i) PVAL_tem]=corr(XS(:,1),Yp);
end

if RHO>0
p_RHO=length(find(RHO<=RHO_tem))/rep;
else
 p_RHO=length(find(RHO>=RHO_tem))/rep;   
end
%permutation test

hold on;plot(XS(:,1),Yhat,'-r');title(['R=',num2str(RHO),'--P=',num2str(p_RHO)])
xlabel('PLS1-score');ylabel('T-value');
print(gcf,[Str,'_gene_image_assocaition.dng'],'-dpng');
savefig(gcf,[Str,'_gene_image_assocoaition.fig']);
close all;
% 

if length(P_gene)>0;
[RHO_P1 PVAL_P1]=corr(X1(:,ind_P(1)),Y);
[RHO_Pend PVAL_Pend]=corr(X1(:,ind_P(end)),Y);
end

if length(N_gene)>0
[RHO_N1 PVAL_N1]=corr(X1(:,ind_N(1)),Y);
[RHO_Nend PVAL_Nend]=corr(X1(:,ind_N(end)),Y);
end



% Gene enrichments:
% This code calculates enrichments of the known gene lists (used for Gandal, DISEASES and GAD). The input gene list must be a list of entrez IDs.
cell_data=importdata('...\data\Gene_data\41467_2020_17051_MOESM8_ESM.xlsx');
raw_cell=cell_data.textdata(2:end,:);
celltype_name={'Astro' 'Endo' 'Micro' 'Neuro' 'Neuro-Ex' 'Neuro-In' 'Oligo' 'OPC'};%Per
clear list_celltype;
list_celltype=cell(length(celltype_name),1);
for i=1:length(celltype_name)     
    if sum((strcmp(celltype_name{i},raw_cell(:,4))))>0
        Index_cell=find(strcmp(celltype_name{i},raw_cell(:,4)));
        for z=1:length(Index_cell)
        list_celltype{i}=[list_celltype{i} raw_cell(Index_cell(z),5:end)];
        
        end
       list_celltype{i}(cellfun(@isempty,list_celltype{i}))=[];
       tem_knowngenes=cellfun(@num2str,list_celltype{i},'un',0);
       list_celltype{i}=unique(tem_knowngenes);
       %[~,k] = unique(cellfun(@char,cellfun(@getByteStreamFromArray,list_celltype{i},'un',0),'un',0));
    end
end

bground=[list_celltype{1} list_celltype{2} list_celltype{3} list_celltype{4} list_celltype{5} list_celltype{6} list_celltype{7} list_celltype{8}];
bground=unique(bground);
%for PLS-
for i=1:length(celltype_name)
    i
knowngenes=list_celltype{i}; % List of known genes (entrez IDs)
knowngenes=unique(knowngenes);
names=Genename(Label<0); % for PLS-genelist 

tem_p=hygepdf(length(intersect(knowngenes,names)):length(intersect(bground,names)),length(bground),length(knowngenes),length(intersect(bground,names)));

Neg_p_gene_overlap(i)=sum(tem_p);
%PLS+
knowngenes=list_celltype{i}; % List of known genes (entrez IDs)
knowngenes=unique(knowngenes);
names=Genename(Label>0); % for PLS-genelist 

tem_p=hygepdf(length(intersect(knowngenes,names)):length(intersect(bground,names)),length(bground),length(knowngenes),length(intersect(bground,names)));

Pos_p_gene_overlap(i)=sum(tem_p);
end

Neg_FDR_p_gene_overlap=mafdr(Neg_p_gene_overlap,'BHFDR','true');
Pos_FDR_p_gene_overlap=mafdr(Pos_p_gene_overlap,'BHFDR','true');

save([Str '.mat'])

Pos_FDR_p_gene_overlap=Pos_FDR_p_gene_overlap';
Neg_FDR_p_gene_overlap=Neg_FDR_p_gene_overlap';



for i=1:length(celltype_name)
    i
knowngenes=list_celltype{i}; % List of known genes (entrez IDs)
knowngenes=unique(knowngenes);
names=Genename(Label<0); % for PLS-genelist 

%tem_p=hygepdf(length(intersect(knowngenes,names)):length(intersect(bground,names)),length(bground),length(knowngenes),length(intersect(bground,names)));

%Neg_p_gene_overlap(i)=sum(tem_p);

Cell_Results_Neg{i,1}=celltype_name{i};
Cell_Results_Neg{i,2}=length(knowngenes);
Cell_Results_Neg{i,3}=length(names);
Cell_Results_Neg{i,4}=length(intersect(knowngenes,names));
Cell_Results_Neg{i,5}=Neg_p_gene_overlap(i);
Cell_Results_Neg{i,6}=Neg_FDR_p_gene_overlap(i);
Cell_Results_Neg{i,7}=Str;

%PLS+
knowngenes=list_celltype{i}; % List of known genes (entrez IDs)
knowngenes=unique(knowngenes);
names=Genename(Label>0); % for PLS-genelist 

%tem_p=hygepdf(length(intersect(knowngenes,names)):length(intersect(bground,names)),length(bground),length(knowngenes),length(intersect(bground,names)));

%Pos_p_gene_overlap(i)=sum(tem_p);
Cell_Results_Pos{i,1}=celltype_name{i};
Cell_Results_Pos{i,2}=length(knowngenes);
Cell_Results_Pos{i,3}=length(names);
Cell_Results_Pos{i,4}=length(intersect(knowngenes,names));
Cell_Results_Pos{i,5}=Pos_p_gene_overlap(i);
Cell_Results_Pos{i,6}=Pos_FDR_p_gene_overlap(i);
Cell_Results_Pos{i,7}=Str;
end

cd(savepath);
xlswrite([Str,'_cell_result.xlsx'],Cell_Results_Neg,'Cell_Results_Neg')
xlswrite([Str,'_cell_result.xlsx'],Cell_Results_Pos,'Cell_Results_Pos')