%-----------------------------------------------------------------------
% Job saved on 10-Sep-2021 12:06:48 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
global T1
cd(T1)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[T1 '\T1W_lesion_filled.nii,1']};    # replace the dtapath for your T1 data in .nii format
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[T1 '\nodif.nii,1']};    #replace the dtapath for your nodif data in .nii format
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
                                                    [T1 '\FIT_ICVF.nii,1']
                                                   [T1 '\FIT_ISOVF.nii,1']
                                                   [T1 '\FIT_OD.nii,1']
                                                   };     #replace the dtapath for your NODDI fitting parameters in .nii format
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'B0_2_T1_';
matlabbatch{2}.spm.spatial.normalise.write.subj.def = {[T1 '\y_T1W_lesion_filled.nii']};
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
                                                          90 90 108];
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 2;
