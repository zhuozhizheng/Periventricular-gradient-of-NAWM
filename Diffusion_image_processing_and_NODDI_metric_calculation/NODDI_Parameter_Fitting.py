import amico
import os
filepath='.../data/NODDI' # the datapath to your preprcessed diffusion data (many subjects)
patient_ID=os.listdir(filepath)

for i in patient_ID:
    if not os.path.exists(os.path.join(filepath,i,"AMICO")):
        print("processing "+i)
        os.chdir('.../data')
        print(os.getcwd())
        amico.core.setup()
        ae = amico.Evaluation("NODDI",i)
        amico.util.fsl2scheme("NODDI/"+i+"/bvals.bval","NODDI/"+i+"/rotate_bvecs.bvec")
        ae.load_data(dwi_filename = "DATA.nii.gz", scheme_filename = "bvals.scheme", mask_filename = "nodif_brain_mask.nii.gz", b0_thr = 0)
        ae.set_model("NODDI")
        ae.generate_kernels()
        ae.load_kernels()
        ae.fit()
        ae.save_results()

