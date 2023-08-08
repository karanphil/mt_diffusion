#! /bin/bash
# cd ~/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
cd ~/data/stockage/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=(sub*);
cd ..;

for dir in "${all_dirs[@]}";
    do echo ${dir};
    # mkdir -p corrected_measures/${dir};
    # mkdir -p corrected_measures/${dir}/masks;
    # mkdir -p corrected_measures/${dir}/measures;
    # mkdir -p corrected_measures/${dir}/plots;
    # mkdir -p corrected_measures/${dir}/results_txt;
    mkdir -p corrected_measures_new_nufo/${dir};
    mkdir -p corrected_measures_new_nufo/${dir}/masks;
    mkdir -p corrected_measures_new_nufo/${dir}/measures;
    mkdir -p corrected_measures_new_nufo/${dir}/plots;
    mkdir -p corrected_measures_new_nufo/${dir}/results_txt;
    # scil_image_math.py union bundles/${dir}/rois/CC_*.nii.gz bundles/${dir}/CC_safe.nii.gz -f;
    # python ~/source/mt_diffusion/scripts/correct_metrics_for_orientation.py FODF_metrics/${dir}/${dir}__peaks.nii.gz FODF_metrics/${dir}/${dir}__peak_values.nii.gz DTI_metrics/${dir}/${dir}__dti_fa.nii.gz wm_mask/${dir}/${dir}__wm_mask.nii.gz FODF_metrics/${dir}/${dir}__nufo.nii.gz corrected_measures/${dir}/ --in_mtr ihMT/${dir}/${dir}__MTR_warped.nii.gz --in_ihmtr ihMT/${dir}/${dir}__ihMTR_warped.nii.gz --in_mtsat ihMT/${dir}/${dir}__MTsat_warped.nii.gz --in_ihmtsat ihMT/${dir}/${dir}__ihMTsat_warped.nii.gz --in_e1 DTI_metrics/${dir}/${dir}__dti_evecs_v1.nii.gz --in_roi bundles/${dir}/CC_safe.nii.gz;
    python ~/Research/source/mt_diffusion/scripts/correct_metrics_for_orientation.py FODF_metrics/${dir}/${dir}__peaks.nii.gz FODF_metrics/${dir}/${dir}__peak_values.nii.gz DTI_metrics/${dir}/${dir}__dti_fa.nii.gz wm_mask/${dir}/${dir}__wm_mask.nii.gz recompute_nufo_with_fix_at02/${dir}__nufo.nii.gz corrected_measures_new_nufo/${dir}/ --in_mtr ihMT/${dir}/${dir}__MTR_warped.nii.gz --in_ihmtr ihMT/${dir}/${dir}__ihMTR_warped.nii.gz --in_mtsat ihMT/${dir}/${dir}__MTsat_warped.nii.gz --in_ihmtsat ihMT/${dir}/${dir}__ihMTsat_warped.nii.gz --in_e1 DTI_metrics/${dir}/${dir}__dti_evecs_v1.nii.gz --in_roi bundles/${dir}/CC_safe.nii.gz;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/mtr_corrected_WB.nii.gz corrected_measures_new_nufo/${dir}/measures/mtr_corrected_WB_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/mtr_corrected_CC.nii.gz corrected_measures_new_nufo/${dir}/measures/mtr_corrected_CC_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/mtsat_corrected_WB.nii.gz corrected_measures_new_nufo/${dir}/measures/mtsat_corrected_WB_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/mtsat_corrected_CC.nii.gz corrected_measures_new_nufo/${dir}/measures/mtsat_corrected_CC_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/ihmtr_corrected_WB.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtr_corrected_WB_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/ihmtr_corrected_CC.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtr_corrected_CC_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/ihmtsat_corrected_WB.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtsat_corrected_WB_nlmeans.nii.gz 1 --processes 8 -f;
    scil_run_nlmeans.py corrected_measures_new_nufo/${dir}/measures/ihmtsat_corrected_CC.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtsat_corrected_CC_nlmeans.nii.gz 1 --processes 8 -f;

    scil_image_math.py lower_threshold_eq wm_mask/${dir}/${dir}__wm_mask.nii.gz 0.9 wm_mask/${dir}/${dir}__wm_real_mask.nii.gz -f;

    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/mtr_corrected_WB_nlmeans.nii.gz ihmt/${dir}/${dir}__MTR_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/mtr_difference_WB_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/mtr_difference_WB_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/mtr_difference_WB_nlmeans.nii.gz --data_type float64 -f;
    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/mtr_corrected_CC_nlmeans.nii.gz ihmt/${dir}/${dir}__MTR_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/mtr_difference_CC_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/mtr_difference_CC_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/mtr_difference_CC_nlmeans.nii.gz --data_type float64 -f;

    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/mtsat_corrected_WB_nlmeans.nii.gz ihmt/${dir}/${dir}__MTsat_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/mtsat_difference_WB_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/mtsat_difference_WB_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/mtsat_difference_WB_nlmeans.nii.gz --data_type float64 -f;
    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/mtsat_corrected_CC_nlmeans.nii.gz ihmt/${dir}/${dir}__MTsat_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/mtsat_difference_CC_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/mtsat_difference_CC_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/mtsat_difference_CC_nlmeans.nii.gz --data_type float64 -f;

    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/ihmtr_corrected_WB_nlmeans.nii.gz ihmt/${dir}/${dir}__ihMTR_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtr_difference_WB_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/ihmtr_difference_WB_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtr_difference_WB_nlmeans.nii.gz --data_type float64 -f;
    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/ihmtr_corrected_CC_nlmeans.nii.gz ihmt/${dir}/${dir}__ihMTR_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtr_difference_CC_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/ihmtr_difference_CC_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtr_difference_CC_nlmeans.nii.gz --data_type float64 -f;

    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/ihmtsat_corrected_WB_nlmeans.nii.gz ihmt/${dir}/${dir}__ihMTsat_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtsat_difference_WB_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/ihmtsat_difference_WB_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtsat_difference_WB_nlmeans.nii.gz --data_type float64 -f;
    scil_image_math.py subtraction corrected_measures_new_nufo/${dir}/measures/ihmtsat_corrected_CC_nlmeans.nii.gz ihmt/${dir}/${dir}__ihMTsat_warped.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtsat_difference_CC_nlmeans.nii.gz -f;
    scil_image_math.py multiplication corrected_measures_new_nufo/${dir}/measures/ihmtsat_difference_CC_nlmeans.nii.gz wm_mask/${dir}/${dir}__wm_real_mask.nii.gz corrected_measures_new_nufo/${dir}/measures/ihmtsat_difference_CC_nlmeans.nii.gz --data_type float64 -f;
done;