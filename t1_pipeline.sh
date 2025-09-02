#!/bin/bash

b0=$3;
fa=$2;
t1=$1;
mask=$4;
sub=$5;
type=$6;

processed_data_dir="full_processing/${sub}/${type}";
cd $processed_data_dir;

# mrconvert -strides 1,2,3 $t1 t1.nii.gz;

###############################################################
# T1 processing (can be done in || to diffusion processing)
###############################################################
# bet T1
echo "Bet";
antsBrainExtraction.sh -d 3 -a t1.nii.gz -u 0 -e ../../../t1_template.nii.gz -m ../../../t1_brain_probability_map.nii.gz -o output;
rm -r output;
scil_volume_math.py convert outputBrainExtractionMask.nii.gz t1_bet_mask.nii.gz --data_type uint8 -f;
scil_volume_math.py multiplication t1.nii.gz t1_bet_mask.nii.gz t1_bet.nii.gz --data_type float32 -f;
# denoise T1
echo "Denoise";
# scil_denoising_nlmeans.py t1_bet.nii.gz t1_bet_nlm.nii.gz --mask_denoise t1_bet_mask.nii.gz --processes 6 --number_coils 1 --basic_sigma -f;
# With the singularity, use this old version of the script:
scil_denoising_nlmeans.py t1_bet.nii.gz t1_bet_nlm.nii.gz 1 --mask t1_bet_mask.nii.gz --processes 6 -f;
# Segment T1
echo "Segmentation";
fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -g -o t1_bet_nlm.nii.gz t1_bet_nlm.nii.gz;

# rename output of T1 segmentation
echo "Rename";
mv t1_bet_nlm_seg_2.nii.gz wm_mask.nii.gz;
mv t1_bet_nlm_seg_1.nii.gz gm_mask.nii.gz;
mv t1_bet_nlm_seg_0.nii.gz csf_mask.nii.gz;
mv t1_bet_nlm_pve_2.nii.gz wm_map.nii.gz;
mv t1_bet_nlm_pve_1.nii.gz gm_map.nii.gz;
mv t1_bet_nlm_pve_0.nii.gz csf_map.nii.gz;
rm -rf t1_bet_nlm_*;

mkdir register_natif; cd register_natif;

echo "Registration";
t1="../t1_bet_nlm.nii.gz";
antsRegistration --dimensionality 3 --float 0\
		 --output [output,outputWarped.nii.gz,outputInverseWarped.nii.gz]\
		 --interpolation Linear --use-histogram-matching 0\
		 --winsorize-image-intensities [0.005,0.995]\
		 --initial-moving-transform [$b0,$t1,1]\
		 --transform Rigid['0.2']\
		 --metric MI[$b0,$t1,1,32,Regular,0.25]\
		 --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
		 --smoothing-sigmas 3x2x1x0\
		 --transform Affine['0.1']\
		 --metric MI[$b0,$t1,1,32,Regular,0.25]\
		 --convergence [500x250x125x100,1e-6,10] --shrink-factors 8x4x2x1\
		 --smoothing-sigmas 3x2x1x0\
		 --transform SyN[0.01,3,0]\
		 --metric MI[$b0,$t1,1,32]\
		 --metric CC[$fa,$t1,1,4]\
		 --convergence [50x25x10,1e-6,10] --shrink-factors 4x2x1\
		 --smoothing-sigmas 3x2x1\
		 --verbose > log_ants.txt
cd ../;

antsApplyTransforms -d 3 -i wm_mask.nii.gz -r $fa -o register_natif/wm_mask.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz;
antsApplyTransforms -d 3 -i gm_mask.nii.gz -r $fa -o register_natif/gm_mask.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz;
antsApplyTransforms -d 3 -i csf_mask.nii.gz -r $fa -o register_natif/csf_mask.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz;
antsApplyTransforms -d 3 -i wm_map.nii.gz -r $fa -o register_natif/wm_map.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz;
antsApplyTransforms -d 3 -i gm_map.nii.gz -r $fa -o register_natif/gm_map.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz;
antsApplyTransforms -d 3 -i csf_map.nii.gz -r $fa -o register_natif/csf_map.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz;

cd register_natif;
scil_volume_math.py multiplication wm_mask.nii.gz $mask wm_mask.nii.gz --data_type uint16 -f;
scil_volume_math.py multiplication gm_mask.nii.gz $mask gm_mask.nii.gz --data_type uint16 -f;
scil_volume_math.py multiplication csf_mask.nii.gz $mask csf_mask.nii.gz --data_type uint16 -f;
scil_volume_math.py convert wm_mask.nii.gz wm_mask.nii.gz --data_type uint8 -f;
scil_volume_math.py convert gm_mask.nii.gz gm_mask.nii.gz --data_type uint8 -f;
scil_volume_math.py convert csf_mask.nii.gz csf_mask.nii.gz --data_type uint8 -f;

cd ../;

rm -r tmp*;