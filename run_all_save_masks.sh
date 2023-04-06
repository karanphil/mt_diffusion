#! /bin/bash
cd /home/pkaran/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=$(ls);

cd ../output_ratios;

for dir in $all_dirs;
    do python ~/source/mt_diffusion/save_angle_bins_mask.py ../DTI_metrics/${dir}/${dir}__dti_evecs_v1.nii.gz ../wm_mask/${dir}/${dir}__wm_mask.nii.gz;
    echo $dir;
done;
