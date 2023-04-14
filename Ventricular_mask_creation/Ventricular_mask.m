%Step1 create CSF group mask
clear all;clc;
sum_data=zeros(121,145,121);

datapath='...\data\NODDI';%path to data
savepath='...\results'%path to save results
cd(datapath);
f=dir;f1=f(3:end);


for i=1:length(f1)
    i
    path=[datapath filesep f1(i).name];
    cd(path);
    tem=spm_vol('wc3_T1W_lesion_filled.nii');
    data=spm_read_vols(tem);
    sum_data=sum_data+data;    
end

sum_data(sum_data<length(f1)*0.01)=0;%remove the voxel that in less than 1% subjects
cd(savepath);
tem.fname='CSF_group_mask.nii';
spm_write_vol(tem,sum_data/length(f1));

%Step2 manually labelled ventricular group mask using ITKSnap and modified

%to remove the potential small isolated clusters (voxel<100)
tem=spm_vol('Ventricular_group_level_1.nii');
data=spm_read_vols(tem);
data(data<0.01)=0;% to confirm the backgroud as zero
[L num]=bwlabeln(data);

for j=1:num
    if length(find(L==j))<100
        L(L==j)=0;
    end
end
L(L>0)=1;
tem.fname='New_Ventricular_group_level_final.nii';
spm_write_vol(tem,L);


%step3 calculate ventricular individual mask
clear all;clc;

datapath='...\data\NODDI';

cd(datapath);
f=dir;f1=f(3:end);

tem=spm_vol('...\results\Ventricular_group_level_final.nii');
mask=spm_read_vols(tem);

for i=1:length(f1)
    i
    path=[datapath filesep f1(i).name];
    cd(path);
    tem=spm_vol('wc3_T1W_lesion_filled.nii);
    data=spm_read_vols(tem);
    data(data<0.5)=0;
    data=data.*mask;
    data(data>0)=1;
    
    [L num]=bwlabeln(data,26);
    
    for j=1:num
        if length(find(L==j))<100
            L(L==j)=0;
        end
    end
    L(L>0)=1;
    
    
    tem.fname='Individual_ven_mask.nii';
    spm_write_vol(tem,L);
end





