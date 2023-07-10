#! /bin/bash
cd ~/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=(sub*);
cd ..;

for dir in "${all_dirs[@]}";
    do echo ${dir};
    scil_run_nlmeans.py corrected_measures/${dir}/measures/mtr_corrected_WB.nii.gz corrected_measures/${dir}/measures/mtr_corrected_WB.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures/${dir}/measures/mtsat_corrected_WB.nii.gz corrected_measures/${dir}/measures/mtsat_corrected_WB.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures/${dir}/measures/ihmtr_corrected_WB.nii.gz corrected_measures/${dir}/measures/ihmtr_corrected_WB.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures/${dir}/measures/ihmtsat_corrected_WB.nii.gz corrected_measures/${dir}/measures/ihmtsat_corrected_WB.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures/${dir}/measures/mtr_corrected_CC.nii.gz corrected_measures/${dir}/measures/mtr_corrected_CC.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures/${dir}/measures/ihmtr_corrected_CC.nii.gz corrected_measures/${dir}/measures/ihmtr_corrected_CC.nii.gz 1 --processes 8 -f;
    done;