% note: this code is used to calculate the periventricular distance map and rings (3mm thick per ring) within the WM (wc2 data with threshold >0.9) 
clear all;clc;
datapath='..data/NODDI';
cd(datapath);f=dir;f1=f(3:end);

for z=1:length(f1)
 
voxel_size=1.5;
path=[datapath filesep f1(z).name];

cd(path);tem=spm_vol('Individual_ven_mask.nii'); %load the individual ventricle mask
CSF0=spm_read_vols(tem);CSF0(CSF0>0)=1;CSF0(CSF0<=0)=0;
CSF=CSF0;


Idx1=find(CSF>0);
[v11,v12,v13]=ind2sub(size(CSF),Idx1);
V1=[v11 v12 v13];

cd(path);tem=spm_vol('wc2_T1W_lesion_filled.nii'); %load the WM segmentation using CAT12
WM=spm_read_vols(tem);WM(WM<0.9)=0;  %threshold the WM segmentation with cutoff of 0,9 to restricted the calculation within WM.

Idx0=find(WM>0);
[v01,v02,v03]=ind2sub(size(WM),Idx0);
V0=[v01 v02 v03];

D=pdist2(V0,V1,'euclidean');% calculate the euclidean distance of beteen the voxels within ventricle and voxels within WM

m=size(V0,1);
tem_dis_map=min(D,[],2);% determine the minmum euclidean distance of beteen the paired voxels within ventricle and WM
Dis_map=zeros(size(WM));


for i=1:m
    
    Dis_map(V0(i,1),V0(i,2),V0(i,3))=tem_dis_map(i);
    
end

Dis_map=Dis_map;
Ring_dis_map=Dis_map;
Max_value=max(max(max(Dis_map)));
Inter_value=2; % set the ring thickness as 3mm corresonding to 2 volex size (1.5mm)
seg_num=round(Max_value/Inter_value);

for i=1:seg_num;
    
    tem_dis_map1=zeros(size(Dis_map));
    tem_dis_map2=zeros(size(Dis_map)); 
    tem_dis_map3=zeros(size(Dis_map));
    
    tem_dis_map1(Dis_map<=i*Inter_value)=1;
    tem_dis_map2(Dis_map>(i-1)*Inter_value)=1;
    
    tem_dis_map3=tem_dis_map1.*tem_dis_map2;
    Ring_dis_map(tem_dis_map3>0)=i;
end
cd(path);
tem=spm_vol('wc2_T1W_lesion_filled.nii');
tem.fname='Distance_map.nii'; %save the raw euclidean distance maps 
spm_write_vol(tem,Dis_map);

tem.fname='Ring_distance_map.nii'; %save the ring index maps 
spm_write_vol(tem,Ring_dis_map);

end