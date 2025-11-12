#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh and
# preprocessing_t1_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------FODF PROCESSING-------------------------
b0_thr_extract_b0=10;

echo "FODF PROCESSING";
for sub in $subs; 
    do echo $sub;

    dwi_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.nii.gz";
    bval_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.bval";
    bvec_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.bvec";
    dwi_mt_on="${target_dir}/${sub}/dwi/dwi_mt_on.nii.gz";
    bval_mt_on="${target_dir}/${sub}/dwi/dwi_mt_on.bval";
    bvec_mt_on="${target_dir}/${sub}/dwi/dwi_mt_on.bvec";
    mask="${target_dir}/${sub}/dwi/b0_brain_mask.nii.gz";
    b0="${target_dir}/${sub}/dwi/b0.nii.gz";
    wm_mask="${target_dir}/${sub}/preprocessing_t1/register_natif/wm_mask.nii.gz";
    csf_mask="${target_dir}/${sub}/preprocessing_t1/register_natif/csf_mask.nii.gz";

	# ----------------------DTI Computation---------------------
    cd ${target_dir}/${sub};
    mkdir -p dti;
    cd ${target_dir}/${sub}/dti;
    # This should already be done in the t1 preprocessing, but just in case.
	if [ ! -f "fa.nii.gz" ]; then
		echo "DTI computation";
		scil_dti_metrics.py $dwi_mt_off $bval_mt_off $bvec_mt_off --mask $mask --not_all --fa fa.nii.gz --md md.nii.gz --rgb rgb.nii.gz -f;
	fi
	fa="${target_dir}/${sub}/dti/fa.nii.gz";
    md="${target_dir}/${sub}/dti/md.nii.gz";

    # ---------------------CSD Computation---------------------
    cd ${target_dir}/${sub};
    mkdir -p fodf;
    cd ${target_dir}/${sub}/fodf;
    if [ ! -f "fodf_mt_off.nii.gz" ]; then
        echo "CSD Computation";
        scil_frf_ssst.py $dwi_mt_off $bval_mt_off $bvec_mt_off frf.txt --mask $mask --mask_wm $wm_mask --roi_radii 15 15 10 -f;
        scil_fodf_ssst.py $dwi_mt_off $bval_mt_off $bvec_mt_off frf.txt fodf_mt_off.nii.gz --mask $mask --processes 8 --sh_order 6 -f;
        scil_fodf_ssst.py $dwi_mt_on $bval_mt_on $bvec_mt_on frf.txt fodf_mt_on.nii.gz --mask $mask --processes 8 --sh_order 6 -f;
    fi
    fodf_mt_off="${target_dir}/${sub}/fodf/fodf_mt_off.nii.gz";
    fodf_mt_on="${target_dir}/${sub}/fodf/fodf_mt_on.nii.gz";

    # ---------------------FODF Metrics Computation---------------------
    # MT-off
    cd ${target_dir}/${sub};
    mkdir fodf_metrics_mt_off;
    cd ${target_dir}/${sub}/fodf_metrics_mt_off;
    if [ ! -f "peaks.nii.gz" ]; then
        echo "FODF metrics on MT-off";
        scil_fodf_max_in_ventricles.py $fodf_mt_off $fa $md --md_threshold 0.0025 --max_value_output max_fodf_in_ventricles.txt --in_mask $csf_mask --use_median -f;
        max_value=$(cat max_fodf_in_ventricles.txt);    
        a_threshold=$(echo 2*${max_value}|bc);
        scil_fodf_metrics.py $fodf_mt_off --mask $mask --abs_peaks_and_values --at $a_threshold -f --processes 8;
    fi
    # MT-on
    cd ${target_dir}/${sub};
    mkdir fodf_metrics_mt_on;
    cd ${target_dir}/${sub}/fodf_metrics_mt_on;
    if [ ! -f "peaks.nii.gz" ]; then
        echo "FODF metrics on MT-on";
        scil_fodf_max_in_ventricles.py $fodf_mt_on $fa $md --md_threshold 0.0025 --max_value_output max_fodf_in_ventricles.txt --in_mask $csf_mask --use_median -f;
        max_value=$(cat max_fodf_in_ventricles.txt);
        a_threshold=$(echo 2*${max_value}|bc);
        scil_fodf_metrics.py $fodf_mt_on --mask $mask --abs_peaks_and_values --at $a_threshold -f --processes 8;
    fi

done;

