#! /bin/bash
cd /home/pkaran/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=$(ls);

cd ../output;

for dir in $all_dirs;
    do python ~/source/mt_diffusion/average_mt_ihmt_by_angle.py ../DTI_metrics/${dir}/${dir}__dti_evecs_v1.nii.gz ../DTI_metrics/${dir}/${dir}__dti_fa.nii.gz ../ihMT/${dir}/${dir}__MTR_warped.nii.gz ../ihMT/${dir}/${dir}__ihMTR_warped.nii.gz ../wm_mask/${dir}/${dir}__wm_mask.nii.gz ../FODF_metrics/${dir}/${dir}__fodf_nufo.nii.gz;
    echo $dir;
done;
