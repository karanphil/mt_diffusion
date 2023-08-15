#! /bin/bash
cd ~/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=(sub*);
cd ..;

for dir in "${all_dirs[@]}";
    do echo ${dir};
    mkdir ~/Samsung/data/MT_Diffusion/myelo_inferno/FODF_metrics/${dir}/new_peaks;
    cd ~/Samsung/data/MT_Diffusion/myelo_inferno/FODF_metrics/${dir}/new_peaks;
    scil_compute_fodf_metrics.py ~/Samsung/data/MT_Diffusion/myelo_inferno/FODF_metrics/${dir}/${dir}__fodf.nii.gz --abs_peaks_and_values --at 0.2 --mask ~/Samsung/data/MT_Diffusion/myelo_inferno/b0_mask/${dir}/${dir}__b0_mask.nii.gz --processes 8 -f;
done;