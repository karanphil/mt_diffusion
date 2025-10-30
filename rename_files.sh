#!/bin/bash

# The first argument of the script is the directory where the original data is stored (full path)
# The second argument of the script is the target directory (full path)
original_dir=$1; # ex: /data/karp2601/stockage/mt-diff-mcgill/for-philippe/mt-diff-10peeps
target_dir=$2; # ex: /home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing
cd $original_dir;
subs=$(ls -d hc*);

for sub in $subs; 
    do echo $sub;
    mkdir -p ${target_dir}/${sub}/renamed_data;
    cp ${sub}/dicom/dicom_64_MPRAGE*6.nii ${target_dir}/${sub}/renamed_data/t1.nii;
    gzip ${target_dir}/${sub}/renamed_data/t1.nii;
    cp ${sub}/dicom/dicom_irl_ep2d_mt_diff*.bval ${target_dir}/${sub}/renamed_data/mt_on_dwi.bval;
    cp ${sub}/dicom/dicom_irl_ep2d_mt_diff*.bvec ${target_dir}/${sub}/renamed_data/mt_on_dwi.bvec;
    cp ${sub}/dicom/dicom_irl_ep2d_mt_diff*.nii ${target_dir}/${sub}/renamed_data/mt_on_dwi.nii;
    gzip ${target_dir}/${sub}/renamed_data/mt_on_dwi.nii;
    cp ${sub}/dicom/dicom_irl_ep2d_nomt_diff_PL_2*.bval ${target_dir}/${sub}/renamed_data/mt_off_dwi.bval;
    cp ${sub}/dicom/dicom_irl_ep2d_nomt_diff_PL_2*.bvec ${target_dir}/${sub}/renamed_data/mt_off_dwi.bvec;
    cp ${sub}/dicom/dicom_irl_ep2d_nomt_diff_PL_2*.nii ${target_dir}/${sub}/renamed_data/mt_off_dwi.nii;
    gzip ${target_dir}/${sub}/renamed_data/mt_off_dwi.nii;
    cp ${sub}/dicom/dicom_irl_ep2d_nomt_diff_PL_b0PA*.nii ${target_dir}/${sub}/renamed_data/mt_off_revb0.nii;
    gzip ${target_dir}/${sub}/renamed_data/mt_off_revb0.nii;

done;