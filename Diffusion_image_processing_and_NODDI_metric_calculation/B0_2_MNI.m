% List of open inputs
clear all;clc
global T1
T1path='...\data\NODDI'; # the datapath to you NODDI data 
cd(T1path);f=dir;f1=f(3:end);

for i=1:length(f1)
    T1=[T1path filesep f1(i).name];
    
    cd(T1); 
    ICVF_file_gz=dir('FIT_ICVF.nii.gz');
    ICVF_file_nii=dir('FIT_ICVF.nii');
    
    if length(ICVF_file_gz)>0&&length(ICVF_file_nii)<1
        FIT_gz=dir('FIT*.nii.gz')
        for j=1:length(FIT_gz)
        gunzip(FIT_gz(j).name);
        end
        
        data_gz=dir('data*.nii.gz');
        for j=1:length(data_gz)
        gunzip(data_gz(j).name);
        end     
    end
    
    y_T1_file=dir('y_T1W_lesion_filled.nii');
    ICVF_file=dir('FIT_ICVF.nii');
    B0_file=dir('nodif.nii');
    
    if length(y_T1_file)==0||length(ICVF_file)==0||length(B0_file)==0
        list{i,1}=f1(i).name;
        list{i,2}=length(y_T1_file);
        list{i,3}=length(ICVF_file);
        list{i,4}=length(B0_file);
        
        continue;
    else
        
    nrun = 1; % enter the number of runs here
    jobfile = {'..../B0_2_MNI_job.m'}; #replace the path of your B0_2_MNI_job path
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    for crun = 1:nrun
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
    end
    
end
