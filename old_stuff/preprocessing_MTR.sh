#!/bin/bash
# À rouler dans /home/local/USHERBROOKE/karp2601/Samsung/data/mt-diff-mcgill
# use : singularity exec -B /media/karp2601/T7/data/mt-diff-mcgill/ ~/Research/containers/scilus_2.0.2_from_docker.sif bash code/preprocessing_MTR.sh

home="/media/karp2601/T7/data";
# home="/home/pkaran/Samsung/data";

processed_data_dir="processed_data/preprocessing_MTR";
mask=$home/mt-diff-mcgill/processed_data/preprocessing_t1/t1_bet_mask.nii.gz;
fa=$home/mt-diff-mcgill/processed_data/dti/fa.nii.gz;
affine=$home/mt-diff-mcgill/processed_data/preprocessing_t1/register_natif/output0GenericAffine.mat;
warp=$home/mt-diff-mcgill/processed_data/preprocessing_t1/register_natif/output1Warp.nii.gz;
cd $processed_data_dir;

# bet 
echo "Bet";
scil_volume_reshape_to_reference.py $home/mt-diff-mcgill/original_data/processed/MTR.nii $mask  MTR.nii.gz -f;
scil_volume_math.py multiplication MTR.nii.gz $mask MTR_bet.nii.gz --data_type float32 -f;
# denoise T1
echo "Denoise";
# scil_denoising_nlmeans.py t1_bet.nii.gz t1_bet_nlm.nii.gz --mask_denoise t1_bet_mask.nii.gz --processes 6 --number_coils 1 --basic_sigma -f;
# With the singularity, use this old version of the script:
scil_denoising_nlmeans.py MTR_bet.nii.gz MTR_bet_nlm.nii.gz 1 --mask $mask --processes 6 -f;

mkdir register_natif;

antsApplyTransforms -d 3 -i MTR_bet_nlm.nii.gz -r $fa -o register_natif/MTR_bet_nlm.nii.gz -n NearestNeighbor -t $affine $warp;
