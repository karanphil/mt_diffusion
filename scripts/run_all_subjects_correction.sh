#! /bin/bash
cd ~/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=(sub*);
cd ..;

for dir in "${all_dirs[@]}";
    do echo ${dir};
    mkdir -p corrected_measures/${dir};
    mkdir -p corrected_measures/${dir}/masks;
    mkdir -p corrected_measures/${dir}/measures;
    mkdir -p corrected_measures/${dir}/plots;
    mkdir -p corrected_measures/${dir}/results_txt;
    scil_image_math.py union bundles/${dir}/rois/CC_*.nii.gz bundles/${dir}/CC_safe.nii.gz -f;
    python ~/Research/source/mt_diffusion/scripts/correct_metrics_for_orientation.py FODF_metrics/${dir}/${dir}__peaks.nii.gz FODF_metrics/${dir}/${dir}__peak_values.nii.gz DTI_metrics/${dir}/${dir}__dti_fa.nii.gz wm_mask/${dir}/${dir}__wm_mask.nii.gz FODF_metrics/${dir}/${dir}__nufo.nii.gz corrected_measures/${dir}/ --in_mtr ihMT/${dir}/${dir}__MTR_warped.nii.gz --in_ihmtr ihMT/${dir}/${dir}__ihMTR_warped.nii.gz --in_mtsat ihMT/${dir}/${dir}__MTsat_warped.nii.gz --in_ihmtsat ihMT/${dir}/${dir}__ihMTsat_warped.nii.gz --in_e1 DTI_metrics/${dir}/${dir}__dti_evecs_v1.nii.gz --in_roi bundles/${dir}/CC_safe.nii.gz;
done;