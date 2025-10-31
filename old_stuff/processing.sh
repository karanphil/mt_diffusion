#!/bin/bash
# À rouler dans /home/local/USHERBROOKE/karp2601/Samsung/data/mt-diff-mcgill
home="/home/local/USHERBROOKE/karp2601";
# home="/home/pkaran";

dwi="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_MTR.nii.gz";
dwi_mt_off="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_mt_off_upsample.nii.gz";
bval="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_mt_off.bval";
bvec="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_mt_off.bvec";
mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/b0_mt_off_brain_mask.nii.gz";
#wm_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1/register_natif/wm_mask.nii.gz";
wm_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1/register_natif/wm_mask.nii.gz";
gm_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1/register_natif/gm_mask.nii.gz";
csf_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1/register_natif/csf_mask.nii.gz";
fa="$home/Samsung/data/mt-diff-mcgill/processed_data/dti/fa.nii.gz";
md="$home/Samsung/data/mt-diff-mcgill/processed_data/dti/md.nii.gz";

wm_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1_mediumvox/register_natif/wm_mask.nii.gz";
gm_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1_mediumvox/register_natif/gm_mask.nii.gz";
csf_mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_t1_mediumvox/register_natif/csf_mask.nii.gz";
fa="$home/Samsung/data/mt-diff-mcgill/processed_data/dti_mediumvox/fa.nii.gz";
md="$home/Samsung/data/mt-diff-mcgill/processed_data/dti_mediumvox/md.nii.gz";

# DTI
# To run before preprocessing T1.
cd processed_data/dti;
scil_dti_metrics.py $dwi_mt_off $bval $bvec --mask $mask -f;
cd ..;

# QBALL
cd processed_data/qball;
scil_qball_metrics.py $dwi $bval $bvec --mask $wm_mask --processes 8 --use_qball --sh_order 6 -f;
cd ..;

# MODF
cd processed_data/modf;
dwi="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_MTR_strides.nii.gz";
mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/b0_mt_off_brain_mask_strides.nii.gz";
bvec="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_mt_off_strides.bvec";
# For 1mm iso, put roi_radii to 30 30 20, 2mm iso is 15 15 10
scil_frf_ssst.py $dwi $bval $bvec frf.txt --mask $mask --mask_wm $wm_mask --roi_radii 30 30 20 -f;

scil_fodf_ssst.py $dwi $bval $bvec frf.txt fodf.nii.gz --mask $mask --processes 8 --sh_order 6 -f;

# Filter mODF
scil_sh_to_rish.py fodf.nii.gz rish -f;
scil_volume_math.py upper_threshold rish0.nii.gz 1 normal_rish0_mask.nii.gz -f;
scil_volume_math.py upper_threshold rish2.nii.gz 1 normal_rish2_mask.nii.gz -f;
scil_volume_math.py upper_threshold rish4.nii.gz 1 normal_rish4_mask.nii.gz -f;
scil_volume_math.py upper_threshold rish6.nii.gz 1 normal_rish6_mask.nii.gz -f;
# scil_volume_math.py multiplication normal_rish0_mask.nii.gz normal_rish2_mask.nii.gz normal_rish4_mask.nii.gz normal_rish6_mask.nii.gz normal_rish_mask.nii.gz -f;
scil_volume_math.py multiplication normal_rish2_mask.nii.gz normal_rish4_mask.nii.gz normal_rish6_mask.nii.gz normal_rish_mask.nii.gz -f;
mrcalc fodf.nii.gz normal_rish_mask.nii.gz -mult fodf_cleaned.nii.gz -force;

# Ceci march aussi pour filtrer, mais moins bien
# scil_volume_math.py upper_threshold_eq sf.nii.gz 2 sf_norm_masks.nii.gz -f;
# mrmath sf_norm_masks.nii.gz min sf_norm_mask.nii.gz -axis 3 -force;
# mrcalc fodf.nii.gz sf_norm_mask.nii.gz -mult fodf_cleaned_sf.nii.gz -force;

scil_fodf_max_in_ventricles.py fodf_cleaned.nii.gz $fa $md --max_value_output max_fodf_in_ventricles.txt -f;
max_value=$(cat max_fodf_in_ventricles.txt);
a_threshold=$(echo 2*${max_value}|bc);

scil_fodf_metrics.py fodf_cleaned.nii.gz --mask $mask --abs_peaks_and_values --at $a_threshold -f --processes 8;

# scil_volume_math.py upper_threshold_eq afd_total_sh0.nii.gz 0.0 positive_sh0_mask.nii.gz -f;
# scil_volume_math.py invert positive_sh0_mask.nii.gz positive_sh0_mask.nii.gz -f;
# mrcalc fodf.nii.gz positive_sh0_mask.nii.gz -mult fodf_cleaned.nii.gz -force;
# mrcalc afd_total_sh0.nii.gz positive_sh0_mask.nii.gz -mult afd_total_sh0_positive.nii.gz -force;
# mrcalc peak_values.nii.gz positive_sh0_mask.nii.gz -mult peak_values_positive_sh0.nii.gz -force;

# scil_volume_math.py upper_threshold_eq afd_total_sh0.nii.gz 0.4 high_sh0_mask.nii.gz -f;
# mrcalc fodf_cleaned.nii.gz high_sh0_mask.nii.gz -mult fodf_cleaned_high.nii.gz -force;

# FODF
cd processed_data/fodf;
dwi_mt_off="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_mt_off_strides.nii.gz";
mask="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/b0_mt_off_brain_mask_strides.nii.gz";
bvec="$home/Samsung/data/mt-diff-mcgill/processed_data/preprocessing_dwi/dwi_mt_off_strides.bvec";
# For 1mm iso, put roi_radii to 30 30 20, 2mm iso is 15 15 10
scil_frf_ssst.py $dwi_mt_off $bval $bvec frf.txt --mask $mask --mask_wm $wm_mask --roi_radii 30 30 20 -f;

scil_fodf_ssst.py $dwi_mt_off $bval $bvec frf.txt fodf.nii.gz --mask $mask --processes 8 --sh_order 6 -f; # Even with wm_mask, some fodf explode. I need to find a way to filter these without removing the fodfs in ventricules... or filter all and give 0 to fodf_metrics.

# Filter fODF
scil_sh_to_rish.py fodf.nii.gz rish -f;
scil_volume_math.py upper_threshold rish0.nii.gz 1 normal_rish0_mask.nii.gz -f;
scil_volume_math.py upper_threshold rish2.nii.gz 1 normal_rish2_mask.nii.gz -f;
scil_volume_math.py upper_threshold rish4.nii.gz 1 normal_rish4_mask.nii.gz -f;
scil_volume_math.py upper_threshold rish6.nii.gz 1 normal_rish6_mask.nii.gz -f;
scil_volume_math.py multiplication normal_rish0_mask.nii.gz normal_rish2_mask.nii.gz normal_rish4_mask.nii.gz normal_rish6_mask.nii.gz normal_rish_mask.nii.gz -f;
mrcalc fodf.nii.gz normal_rish_mask.nii.gz -mult fodf_cleaned.nii.gz -force;

scil_fodf_max_in_ventricles.py fodf_cleaned.nii.gz $fa $md --max_value_output max_fodf_in_ventricles.txt -f;
max_value=$(cat max_fodf_in_ventricles.txt);
a_threshold=$(echo 2*${max_value}|bc);
scil_fodf_metrics.py fodf_cleaned.nii.gz --mask $mask --abs_peaks_and_values --at $a_threshold -f --processes 8;

# Tracking
cd processed_data/tracking;
scil_tracking_pft_maps.py $wm_mask $gm_mask $csf_mask --include map_include.nii.gz --exclude map_exclude.nii.gz --interface interface.nii.gz -f;
scil_volume_math.py convert $wm_mask pft_seeding_mask.nii.gz --data_type uint8 -f;
scil_volume_math.py union pft_seeding_mask.nii.gz interface.nii.gz pft_seeding_mask.nii.gz --data_type uint8 -f;
scil_tracking_pft.py ../fodf/fodf.nii.gz pft_seeding_mask.nii.gz map_include.nii.gz map_exclude.nii.gz tmp.trk --algo prob --npv 10 --seed 0 --step 0.5 --theta 20 --min_length 20 --max_length 200 --particles 15 --back 2 --forward 1 --compress 0.2 --sh_basis descoteaux07 -f;
scil_tractogram_remove_invalid.py tmp.trk pft_tracking.trk --remove_single_point -f;

# Step for re-striding (in case of tournier07 basis needed)
# Before performing CSD, all input files must be restrided, and the bvecs flips accordingly
mrconvert -strides 1,2,3,4 DWI DWI_restride;
mrconvert -strides 1,2,3,4 mask mask_restride;
scil_gradients_modify_axes.py dwi_mt_off.bvec dwi_mt_off_strides.bvec -1 2 3; # in case of -1 2 3 4 strides
# CSD, then convert
scil_sh_convert.py fodf.nii.gz fodf_tournier.nii.gz descoteaux07_legacy tournier07 --processes 8 -f;