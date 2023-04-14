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
%

[data2 txt2 raw2]=xlsread('gradient_ring.xlsx',1);
cell_data_ring=importdata('gradient_ring.xlsx');
raw_cell_ring=cell_data_ring.textdata(1,2:end);

genename=raw_cell_ring;
%  for i=2:size(raw2,2)
% genename{1,i-1}=cellstr(raw2{1,i});
%  end

col_start=1;col_end=10;

for disesae=1:5
    switch disease
        case 1
            Y=data1(1,:)';X=data2(col_start:col_end,2:end);Str='HC';
        case 2
            Y=data1(2,:)';X=data2(col_start:col_end,2:end);Str='AD';
        case 3
            Y=data1(3,:)';X=data2(col_start:col_end,2:end);Str='PD';
        case 4
            Y=data1(4,:)';X=data2(col_start:col_end,2:end);Str='SVD';
        case 5
            Y=data1(5,:)';X=data2(col_start:col_end,2:end);Str='MS';
    end
    
    
    Y(isnan(X(:,1))==1,:)=[];X(isnan(X(:,1))==1,:)=[];
    
    
    cd(savepath);
    disp('Gene screening by PLS ....');
    Y(isnan(X(:,1))==1,:)=[];X(isnan(X(:,1))==1,:)=[];
    
    Y0=zscore(Y);X0=zscore(X,0,1);
    
    Per_Y0=perms(Y0);
    
    for per=1:1000   %permutaion times 1000 times considering the ring number of 10
        
        Y0=Per_Y0(per,:)';
        
        [XL YL XS YS BETA PCTVAR MSE stats]=plsregress(X0,Y0);
        
        cd(savepath);
        
        % permutation testing to assess significance of PLS result as a function of
        % the number of components (dim) included:
        rep=5000;dim=1;
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
        Z_W=stats.W(:,1)./std(Per_W,0,2);%Z-score 
        Z_W_5=prctile(Z_W,5,2);
        Z_W_95=prctile(Z_W,95,2);
        
        
        parfor i=1:length(Z_W);
            [H_W P_W(i,1)]=ztest(Z_W(i,1),0,1);%calculation of z and p according to Z-score
        end
        
        [FDR_W Q_W]=mafdr(P_W,'Lambda',[0.01:0.01:0.95]);
        
        Z_weight=stats.W(:,1)./std(XS(:,1));
        
        [Z_W_sort Ind]=sort(Z_W,'ascend');
      
        genename1=genename(Ind);
        Z_weight=Z_weight(Ind);
        X1=X(:,Ind);
        
     
     
        ind_P=find(Z_W_sort>5);% threshold 5
        ind_N=find(Z_W_sort<-5);% threshold 5
        
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
          
        end
        
        
        cd(savepath);xlswrite([Str '_Genename_column' num2str(per) '.csv'],list_gene,1);
        
    end
    
    
end
