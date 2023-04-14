input_path=.../NODDI # the datapath to you NODDI data in nii.gz format (many subjects) 

cd $input_path
for i in *
do
echo $i

cd $i

# note: this additional code is for the NODDI data that with reverse phase-encoding 
#fslroi NODDI_AP.nii.gz AP 0 1
#fslroi NODDI_PA.nii.gz PA 0 1
#fslmerge -t AP_PA AP.nii.gz PA.nii.gz
#topup --imain=AP_PA.nii.gz --datain=acq_param.txt --config=b02b0.cnf --out=AP_PA_topup
#applytopup --imain=NODDI_AP.nii.gz --inindex=1 --datain=acq_param.txt --topup=AP_PA_topup --method=jac --out=AP_cor
#fslroi AP_cor AP_1stVol 0 1
#bet AP_1stVol.nii.gz AP_brain -m -f 0.2
#cp AP_cor.nii.gz data.nii.gz


eddy_correct data.nii.gz DATA.nii.gz 0 
cp *.bvals.bval bvals.bval
cp *.bvecs.bvec bvecs.bvec

dt_rotate_bvecs $i/bvecs.bvec $i/rotate_bvecs.bvec $i/DATA.ecclog

fslroi DATA.nii.gz nodif.nii.gz 0 1

bet nodif.nii.gz nodif_brain.nii.gz -m -f 0.2

cd ..
done
