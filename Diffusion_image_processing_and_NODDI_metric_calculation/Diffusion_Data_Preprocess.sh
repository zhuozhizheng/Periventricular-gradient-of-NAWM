input_path=.../NODDI # the datapath to you NODDI data in nii.gz format (many subjects) 

cd $input_path
for i in *
do
echo $i

cd $i

cp *.bvals.bval bvals.bval
cp *.bvecs.bvec bvecs.bvec

fslroi NODDI_AP.nii.gz AP 0 1
fslroi NODDI_PA.nii.gz PA 0 1
fslmerge -t AP_PA AP.nii.gz PA.nii.gz
topup --imain=AP_PA.nii.gz --datain=acq_param.txt --config=b02b0.cnf --out=AP_PA_topup
applytopup --imain=NODDI_AP.nii.gz --inindex=1 --datain=acq_param.txt --topup=AP_PA_topup --method=jac --out=AP_cor
fslroi AP_cor AP_1stVol 0 1
bet AP_1stVol.nii.gz AP_brain -m -f 0.2
cp AP_cor.nii.gz data.nii.gz

indx=""
for ((i=1; i<=48; i+=1)); do indx="$indx 1"; done
echo $indx > index.txt

eddy --imain=data.nii.gz --mask=AP_brain_mask --acqp=acq_param.txt --index=index.txt --bvecs=bvecs --bvals=bvals --topup=AP_PA_topup --out=eddy_corrected_data

cp eddy_corrected_data.nii.gz DATA.nii.gz

cp eddy_corrected_data.rotated_bvecs bvecs.bvec

cd ..
done

