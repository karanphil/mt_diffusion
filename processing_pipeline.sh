#!/bin/bash

# The script assumes that preprocessing_pipeline.sh has already been run.

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

    dwi_for_tractography="${target_dir}/${sub}/dwi_for_tractography/dwi_for_tractography.nii.gz";
    mask_for_tractography="${target_dir}/${sub}/dwi_for_tractography/b0_for_tractography_brain_mask.nii.gz";
    b0_for_tractography="${target_dir}/${sub}/dwi_for_tractography/b0_for_tractography_brain.nii.gz";