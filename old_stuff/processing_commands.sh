#!/bin/bash
# À rouler dans /home/local/USHERBROOKE/karp2601/Samsung/data/mt-diff-mcgill

# in qball
scil_qball_metrics.py ../dwi_mt_off_norm.nii.gz ../dwi_mt_off.bval ../dwi_mt_off.bvec --mask ../t1/register_natif/wm_mask.nii.gz --processes 8 --not_all --peaks peaks_mt_off.nii.gz --sh sh_mt_off.nii.gz --peak_indices peaks_indices_mt_off.nii.gz --a_power a_power_mt_off.nii.gz --gfa gfa_mt_off.nii.gz --nufo nufo_mt_off.nii.gz -f --use_qball;

# in fodf
scil_frf_ssst.py ../../dwi_mt_off_norm_clip.nii.gz ../../dwi_mt_off.bval ../../dwi_mt_off.bvec frf.txt --mask ../../b0_mt_off_brain_mask.nii.gz --mask_wm ../../t1/register_natif/wm_mask.nii.gz -f;

scil_fodf_ssst.py ../../dwi_mt_off_norm_clip.nii.gz ../../dwi_mt_off.bval ../../dwi_mt_off.bvec frf.txt fodf_mt_off.nii.gz --mask ../../t1/register_natif/wm_mask.nii.gz --processes 8 --sh_order 6 -f;

scil_fodf_max_in_ventricles.py fodf_MTR_WB.nii.gz ../../dti/fa.nii.gz ../../dti/md.nii.gz --max_value_output max_fodf_in_ventricles.txt -f;
max_value=$(cat max_fodf_in_ventricles.txt);
a_threshold=$(echo 2*${max_value}|bc);

scil_fodf_metrics.py ../fodf_MTR.nii.gz --mask ../../../b0_mt_off_brain_mask.nii.gz --abs_peaks_and_values --at $a_threshold -f --processes 8;