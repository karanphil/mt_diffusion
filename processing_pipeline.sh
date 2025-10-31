#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh and
# preprocessing_t1_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------DIFFUSION PRE-PROCESSING-------------------------
b0_thr_extract_b0=10;

echo "DIFFUSION PRE-PROCESSING";
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

    # ----------------------Powder Average MTR Calculation---------------------
    cd ${target_dir}/${sub};
    mkdir -p powder_average;
    cd ${target_dir}/${sub}/powder_average;
    if [ ! -f "powder_averaged_mtr.nii.gz" ]; then
        echo "Powder Average MTR Calculation";
        scil_volume_math subtraction $dwi_mt_off $dwi_mt_on powder_averaged_mtr.nii.gz --data_type float32 -f;
        mrcalc powder_averaged_mtr.nii.gz $dwi_mt_off -div powder_averaged_mtr.nii.gz -force;
        mrcalc powder_averaged_mtr.nii.gz $mask -mult powder_averaged_mtr.nii.gz -force;
        scil_dwi_extract_b0 powder_averaged_mtr.nii.gz $bval_mt_off $bvec_mt_off b0_mtr.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
        scil_volume_math lower_clip b0_mtr.nii.gz 0 b0_mtr.nii.gz -f;
        scil_volume_math upper_clip b0_mtr.nii.gz 1 b0_mtr.nii.gz -f;
        scil_dwi_powder_average powder_averaged_mtr.nii.gz $bval powder_averaged_mtr.nii.gz --mask $mask -f;
        scil_volume_math lower_clip powder_averaged_mtr.nii.gz 0 powder_averaged_mtr.nii.gz -f;
        scil_volume_math upper_clip powder_averaged_mtr.nii.gz 1 powder_averaged_mtr.nii.gz -f;
    fi

	# ----------------------DTI Computation---------------------
    cd ${target_dir}/${sub};
    mkdir -p dti;
    cd ${target_dir}/${sub}/dti;
    # This should already be done in the t1 preprocessing, but just in case.
	if [ ! -f "fa.nii.gz" ]; then
		echo "DTI computation";
		scil_dti_metrics $dwi_mt_off $bval_mt_off $bvec_mt_off --mask $mask --not_all --fa fa.nii.gz --md md.nii.gz --rgb rgb.nii.gz -f;
	fi
	fa="${target_dir}/${sub}/dti/fa.nii.gz";
    md="${target_dir}/${sub}/dti/md.nii.gz";

    

    done;

