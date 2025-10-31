#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh has already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------T1 PRE-PROCESSING----------------------------
echo "DIFFUSION PRE-PROCESSING";
for sub in $subs;
    do echo $sub;
	t1="${target_dir}/${sub}/renamed_data/t1.nii.gz";
	dwi="${target_dir}/${sub}/dwi/dwi_mt_off.nii.gz";
    bval="${target_dir}/${sub}/dwi/dwi_mt_off.bval";
    bvec="${target_dir}/${sub}/dwi/dwi_mt_off.bvec";
    mask="${target_dir}/${sub}/dwi/b0_brain_mask.nii.gz";
	b0="${target_dir}/${sub}/dwi/b0.nii.gz";

	# ----------------------DTI Computation---------------------
    cd ${target_dir}/${sub};
    mkdir -p dti;
    cd ${target_dir}/${sub}/dti;
	if [ ! -f "fa.nii.gz" ]; then
		echo "DTI computation";
			# !!! Some DTI will crash because of NaNs inside the mask. !!!
			# This is then needed in python in the dwi folder:
			# '
			# import nibabel as nib
			# import numpy as np
			# vol = nib.load("dwi_mt_off.nii.gz")
			# dwi = vol.get_fdata()
			# # mask_vol = nib.load("b0_mt_off_brain_mask.nii.gz")
			# # mask = mask_vol.get_fdata().astype(bool)
			# new_dwi = np.where(np.isnan(dwi), 0, dwi)
			# nib.save(nib.Nifti1Image(new_dwi, vol.affine), "dwi_mt_off.nii.gz")
			# '
		scil_dti_metrics $dwi $bval $bvec --mask $mask --not_all --fa fa.nii.gz --md md.nii.gz --rgb rgb.nii.gz -f;
	fi
	fa="${target_dir}/${sub}/dti/fa.nii.gz";

	# ----------------------Re-stride---------------------------
	cd ${target_dir}/${sub};
    mkdir -p preprocessing_t1;
    cd ${target_dir}/${sub}/preprocessing_t1;
	if [ ! -f "t1.nii.gz" ]; then
		echo "Re-stride";
		mrconvert -strides 1,2,3 $t1 t1.nii.gz -force;
	fi

	# ----------------------Bet---------------------------
	if [ ! -f "t1_bet.nii.gz" ]; then
		echo "Bet";
		antsBrainExtraction.sh -d 3 -a t1.nii.gz -u 0 -e ../../../t1_template.nii.gz -m ../../../t1_brain_probability_map.nii.gz -o output;
		scil_volume_math convert outputBrainExtractionMask.nii.gz t1_bet_mask.nii.gz --data_type uint8 -f;
		scil_volume_math multiplication t1.nii.gz t1_bet_mask.nii.gz t1_bet.nii.gz --data_type float32 -f;
	fi

	# ---------------------Denoising---------------------
	if [ ! -f "t1_bet_nlm.nii.gz" ]; then
		echo "Denoise";
		scil_denoising_nlmeans t1_bet.nii.gz t1_bet_nlm.nii.gz --mask_denoise t1_bet_mask.nii.gz --processes 6 --number_coils 1 --basic_sigma -f;
	fi
	t1="${target_dir}/${sub}/preprocessing_t1/t1_bet_nlm.nii.gz";

	# ---------------------Segmentation---------------------
	if [ ! -f "wm_mask.nii.gz" ]; then
		echo "Segmentation";
		fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -g -o t1_bet_nlm.nii.gz t1_bet_nlm.nii.gz;
		mv t1_bet_nlm_seg_2.nii.gz wm_mask.nii.gz;
		mv t1_bet_nlm_seg_1.nii.gz gm_mask.nii.gz;
		mv t1_bet_nlm_seg_0.nii.gz csf_mask.nii.gz;
		mv t1_bet_nlm_pve_2.nii.gz wm_map.nii.gz;
		mv t1_bet_nlm_pve_1.nii.gz gm_map.nii.gz;
		mv t1_bet_nlm_pve_0.nii.gz csf_map.nii.gz;
		rm -rf t1_bet_nlm_*;
	fi

	# ---------------------Registration T1 to DWI---------------------
	mkdir -p register_natif;
	cd register_natif;
	if [ ! -f "output1Warp.nii.gz" ]; then
		echo "Registration";
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
	fi

	# ---------------------Apply Transforms---------------------
	if [ ! -f "wm_mask.nii.gz" ]; then
		echo "Apply transforms";
		antsApplyTransforms -d 3 -i ../wm_mask.nii.gz -r $fa -o wm_mask.nii.gz -n NearestNeighbor -t output0GenericAffine.mat output1Warp.nii.gz;
		antsApplyTransforms -d 3 -i ../gm_mask.nii.gz -r $fa -o gm_mask.nii.gz -n NearestNeighbor -t output0GenericAffine.mat output1Warp.nii.gz;
		antsApplyTransforms -d 3 -i ../csf_mask.nii.gz -r $fa -o csf_mask.nii.gz -n NearestNeighbor -t output0GenericAffine.mat output1Warp.nii.gz;
		antsApplyTransforms -d 3 -i ../wm_map.nii.gz -r $fa -o wm_map.nii.gz -n NearestNeighbor -t output0GenericAffine.mat output1Warp.nii.gz;
		antsApplyTransforms -d 3 -i ../gm_map.nii.gz -r $fa -o gm_map.nii.gz -n NearestNeighbor -t output0GenericAffine.mat output1Warp.nii.gz;
		antsApplyTransforms -d 3 -i ../csf_map.nii.gz -r $fa -o csf_map.nii.gz -n NearestNeighbor -t output0GenericAffine.mat output1Warp.nii.gz;
		scil_volume_math multiplication wm_mask.nii.gz $mask wm_mask.nii.gz --data_type uint16 -f;
		scil_volume_math multiplication gm_mask.nii.gz $mask gm_mask.nii.gz --data_type uint16 -f;
		scil_volume_math multiplication csf_mask.nii.gz $mask csf_mask.nii.gz --data_type uint16 -f;
		scil_volume_math convert wm_mask.nii.gz wm_mask.nii.gz --data_type uint8 -f;
		scil_volume_math convert gm_mask.nii.gz gm_mask.nii.gz --data_type uint8 -f;
		scil_volume_math convert csf_mask.nii.gz csf_mask.nii.gz --data_type uint8 -f;
	fi

done;
