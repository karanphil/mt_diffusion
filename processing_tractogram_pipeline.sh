#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh and
# preprocessing_t1_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument indicates whether to use gpu or not for tractography
use_gpu=$2; # "true" to use gpu, "false" to use cpu
if [ "$use_gpu" == "true" ]; then
    gpu_option="--use_gpu";
else
    gpu_option="";
fi

cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------Tractogram PROCESSING-------------------------
b0_thr_extract_b0=10;

echo "Tractogram PROCESSING";
for sub in $subs; 
    do echo $sub;

    dwi_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.nii.gz";
    bval_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.bval";
    bvec_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.bvec";
    mask="${target_dir}/${sub}/dwi/b0_brain_mask.nii.gz";
    b0="${target_dir}/${sub}/dwi/b0.nii.gz";
    wm_mask="${target_dir}/${sub}/preprocessing_t1/register_natif/wm_mask.nii.gz";
    csf_mask="${target_dir}/${sub}/preprocessing_t1/register_natif/csf_mask.nii.gz";
    fa="${target_dir}/${sub}/dti/fa.nii.gz";
    md="${target_dir}/${sub}/dti/md.nii.gz";

    # ---------------------Resample-----------------------
    cd ${target_dir}/${sub};
    mkdir -p dwi_for_tractography;
    cd ${target_dir}/${sub}/dwi_for_tractography;
    if [ ! -f "csf_mask.nii.gz" ]; then
        echo "Resample";
        scil_volume_resample.py $dwi_mt_off dwi.nii.gz --voxel_size 1 -f;
        scil_volume_resample.py $b0 b0.nii.gz --voxel_size 1 -f;
        scil_volume_resample.py $fa fa.nii.gz --voxel_size 1 -f;
        scil_volume_resample.py $md md.nii.gz --voxel_size 1 -f;
        scil_volume_resample.py $mask mask.nii.gz --voxel_size 1 --interp nn -f;
        scil_volume_resample.py $wm_mask wm_mask.nii.gz --voxel_size 1 --interp nn -f;
        scil_volume_resample.py $csf_mask csf_mask.nii.gz --voxel_size 1 --interp nn -f;
        cp $bval_mt_off dwi.bval;
        cp $bvec_mt_off dwi.bvec;
    fi
    dwi="${target_dir}/${sub}/dwi_for_tractography/dwi.nii.gz";
    bval="${target_dir}/${sub}/dwi_for_tractography/dwi.bval";
    bvec="${target_dir}/${sub}/dwi_for_tractography/dwi.bvec";
    mask="${target_dir}/${sub}/dwi_for_tractography/mask.nii.gz";
    wm_mask="${target_dir}/${sub}/dwi_for_tractography/wm_mask.nii.gz";
    csf_mask="${target_dir}/${sub}/dwi_for_tractography/csf_mask.nii.gz";
    b0="${target_dir}/${sub}/dwi_for_tractography/b0.nii.gz";
    fa="${target_dir}/${sub}/dwi_for_tractography/fa.nii.gz";
    md="${target_dir}/${sub}/dwi_for_tractography/md.nii.gz";

    # ---------------------CSD Computation---------------------
    cd ${target_dir}/${sub};
    mkdir -p fodf_for_tractography;
    cd ${target_dir}/${sub}/fodf_for_tractography;
    if [ ! -f "fodf.nii.gz" ]; then
        echo "CSD for tractography";
        scil_frf_ssst.py $dwi $bval $bvec frf.txt --mask $mask --mask_wm $wm_mask --roi_radii 30 30 20 -f;
        scil_fodf_ssst.py $dwi $bval $bvec frf.txt fodf.nii.gz --mask $mask --processes 8 --sh_order 6 -f;
    fi
    fodf="${target_dir}/${sub}/fodf_for_tractography/fodf.nii.gz";

    # ---------------------Tractography---------------------
    cd ${target_dir}/${sub};
    mkdir -p tractography;
    cd ${target_dir}/${sub}/tractography;
    if [ ! -f "local_tracking.nii.gz" ]; then
        echo "Tractography";
        scil_tracking_local.py $fodf $wm_mask $wm_mask local_tracking.trk $gpu_option --npv 10 -f;
        scil_tractogram_remove_invalid.py local_tracking.trk local_tracking.trk --remove_single_point -f;
    fi
    tractogram="${target_dir}/${sub}/tractography/local_tracking.trk";

    # ---------------------SIFT2---------------------
    if [ ! -f "sift2_weights.nii.gz" ]; then
        echo "SIFT2";
        fodf_tournier="${target_dir}/${sub}/fodf_tournier.nii.gz";
        scil_sh_convert.py $fodf $fodf_tournier descoteaux07_legacy tournier07 -f;
        tractogram_tck="${target_dir}/${sub}/local_tracking.tck";
        scil_tractogram_convert.py $tractogram $tractogram_tck -f;
        tcksift2 $tractogram_tck $fodf_tournier sift2_weights.txt -force;
        scil_tractogram_dps_math.py $tractogram import "sift2" --in_dps_file sift2_weights.txt --out_tractogram $tractogram -f;
        rm $tractogram_tck $fodf_tournier;
    fi

done;