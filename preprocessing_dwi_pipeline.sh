#!/bin/bash

# The script assumes that all subjects folders have been created previously in the target directory
# and that the renamed data is inside a folder named "renamed_data" in each subject folder,
# as done in the rename_files.sh script (creates all folders and renames the data).

# The first argument of the script is the target directory (full path).
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------DIFFUSION PRE-PROCESSING-------------------------
b0_thr_extract_b0=10;
echo "DIFFUSION PRE-PROCESSING";
for sub in $subs; 
    do echo $sub;
    cd ${target_dir}/${sub};
    mkdir -p preprocessing_dwi;
    cd ${target_dir}/${sub}/preprocessing_dwi;
    renamed_data_dir="../renamed_data";

    dwi_mt_off="${renamed_data_dir}/mt_off_dwi.nii.gz";
    dwi_mt_on="${renamed_data_dir}/mt_on_dwi.nii.gz";
    bval="${renamed_data_dir}/mt_off_dwi.bval";
    bvec="${renamed_data_dir}/mt_off_dwi.bvec";
    rev_b0="${renamed_data_dir}/mt_off_revb0.nii.gz";

    # ---------------------Denoising--------------------
    if [ ! -f "mt_off_dwi_dn.nii.gz" ]; then
        echo "Denoising";
        # Merge mt_off and mt_on for denoising
        mrcat -axis 3 $dwi_mt_off $dwi_mt_on dwi_merge.nii.gz -f;
        dwidenoise dwi_merge.nii.gz dwi_merge_dn.nii.gz -f;
        rm dwi_merge.nii.gz;
        mrconvert -coord 3 0:30 dwi_merge_dn.nii.gz mt_off_dwi_dn.nii.gz -force;
        mrconvert -coord 3 31:end dwi_merge_dn.nii.gz mt_on_dwi_dn.nii.gz -force;
        rm dwi_merge_dn.nii.gz;
    fi
    dwi_mt_off="mt_off_dwi_dn.nii.gz";
    dwi_mt_on="mt_on_dwi_dn.nii.gz";

    # --------------------MT-on extension with b0--------------------
    if [ ! -f "mt_on_dwi_extended.nii.gz" ]; then
        echo "Extract b0";
        scil_dwi_extract_b0 $dwi_mt_off $bval $bvec b0_mt_off.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
        # Add b0 MT-off to MT-on DWI
        scil_volume_math concatenate b0_mt_off.nii.gz $dwi_mt_on mt_on_dwi_extended.nii.gz  -f;
        sed -e 's/^/0 /' $bvec > mt_on_dwi_extended.bvec;
        sed -e 's/^/0 /' $bval > mt_on_dwi_extended.bval;
    fi
    dwi_mt_on_ext="mt_on_dwi_extended.nii.gz";
    bvec_ext="mt_on_dwi_extended.bvec";
    bval_ext="mt_on_dwi_extended.bval";
    b0_mt_off="b0_mt_off.nii.gz";

    # --------------------Topup--------------------
    if [ ! -f "topup_fieldcoef.nii.gz" ]; then
        echo "Topup";
        scil_dwi_prepare_topup_command $b0_mt_off $rev_b0 --out_script --out_prefix topup -f;
        sh topup.sh;
        mrconvert corrected_b0s.nii.gz b0_corrected.nii.gz -coord 3 0 -axes 0,1,2 -nthreads 1 -force;
        bet b0_corrected.nii.gz b0_mt_off_brain -m -R;
    fi
    brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

    # --------------------Eddy--------------------
    if [ ! -f "dwi_mt_off.eddy_rotated_bvecs" ]; then
        echo "Eddy";
        # MT-off
        scil_dwi_prepare_eddy_command $dwi_mt_off $bval $bvec $brain_mask_mt_off --eddy_cmd eddy\
                    --b0_thr $b0_thr_extract_b0\
                    --out_script --fix_seed\
                    --lsr_resampling --slice_drop_correction\
                    --topup topup\
                    --out_prefix dwi_mt_off -f;
        sh eddy.sh;
        cp $bval dwi_mt_off.bval;
        cp dwi_mt_off.eddy_rotated_bvecs dwi_mt_off.bvec;
        # MT-on
        scil_dwi_prepare_eddy_command $dwi_mt_on_ext $bval_ext $bvec_ext $brain_mask_mt_off --eddy_cmd eddy\
                    --b0_thr $b0_thr_extract_b0\
                    --out_script --fix_seed\
                    --lsr_resampling --slice_drop_correction\
                    --topup topup\
                    --out_prefix dwi_mt_on -f;
        sh eddy.sh;
        # Remove the first volume from MT-on (b0 added for eddy and bias correction)
        sed -E 's/^[[:space:]]*[-+]?[0-9]+(\.[0-9]*)?[eE][-+]?[0-9]+[[:space:]]*//' dwi_mt_on.eddy_rotated_bvecs > dwi_mt_on.bvec;
        mrconvert -coord 3 1:31 dwi_mt_on.nii.gz dwi_mt_on.nii.gz -force;
    fi
    dwi_mt_off="dwi_mt_off.nii.gz";
    bval_mt_off="dwi_mt_off.bval";
    bvec_mt_off="dwi_mt_off.bvec";
    dwi_mt_on="dwi_mt_on.nii.gz";
    bval_mt_on="dwi_mt_off.bval";
    bvec_mt_on="dwi_mt_on.bvec";

    # ---------------------Bet-----------------------
    if [ ! -f "b0_mt_off_brain.nii.gz" ]; then
        echo "Extract b0";
        scil_dwi_extract_b0 $dwi_mt_off $bval_mt_off $bvec_mt_off b0_mt_off.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;

        echo "Bet b0";
        bet b0_mt_off.nii.gz b0_mt_off_brain -m;
    fi
    brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

    # ---------------------DWI bias correction---------------------
    if [ ! -f "bias_field_mt_off.nii.gz" ]; then
        echo "Bias correction";
        # MT-off
        dwibiascorrect ants $dwi_mt_off $dwi_mt_off -fslgrad $bvec_mt_off $bval_mt_off -mask $brain_mask_mt_off -bias bias_field_mt_off.nii.gz -force;
        # MT-on, apply the same bias field as MT-off
        scil_dwi_apply_bias_field $dwi_mt_on bias_field_mt_off.nii.gz $dwi_mt_on --mask $brain_mask_mt_off -f;
    fi

    # ---------------------Bet-----------------------
    if [ ! -f "b0_mt_off_brain.nii.gz" ]; then
        echo "Extract b0";
        scil_dwi_extract_b0 $dwi_mt_off $bval_mt_off $bvec_mt_off b0_mt_off.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;

        echo "Bet b0";
        bet b0_mt_off.nii.gz b0_mt_off_brain -m;
    fi
    b0_mt_off="b0_mt_off.nii.gz";
    brain_mt_off="b0_mt_off_brain.nii.gz";
    brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

    # ---------------------Normalization-----------------------
    if [ ! -f "dwi_mt_off_norm.nii.gz" ]; then
        # MT-off
        mrcalc $dwi_mt_off $brain_mt_off -div dwi_mt_off_norm.nii.gz -force;
        mrcalc dwi_mt_off_norm.nii.gz $brain_mask_mt_off -mult dwi_mt_off_norm.nii.gz -force;
        scil_volume_math upper_clip dwi_mt_off_norm.nii.gz 1 dwi_mt_off_norm.nii.gz -f;
        scil_volume_math lower_clip dwi_mt_off_norm.nii.gz 0 dwi_mt_off_norm.nii.gz -f;
        #MT-on
        mrcalc $dwi_mt_on $brain_mt_off -div dwi_mt_on_norm.nii.gz -force;
        mrcalc dwi_mt_on_norm.nii.gz $brain_mask_mt_off -mult dwi_mt_on_norm.nii.gz -force;
        scil_volume_math upper_clip dwi_mt_on_norm.nii.gz 1 dwi_mt_on_norm.nii.gz -f;
        scil_volume_math lower_clip dwi_mt_on_norm.nii.gz 0 dwi_mt_on_norm.nii.gz -f;
    fi
    dwi_mt_off="dwi_mt_off_norm.nii.gz";
    dwi_mt_on="dwi_mt_on_norm.nii.gz";

    # ---------------------Resample-----------------------
    if [ ! -f "dwi_mt_off_upsample.nii.gz" ]; then
        echo "Resample";
        # MT-off
        scil_volume_resample $dwi_mt_off dwi_mt_off_upsample.nii.gz --voxel_size 2 -f;
        # MT-on
        scil_volume_resample $dwi_mt_on dwi_mt_on_upsample.nii.gz --voxel_size 2 -f;
    fi
    dwi_mt_off="${target_dir}/${sub}/preprocessing_dwi/dwi_mt_off_upsample.nii.gz";
    bval_mt_off="${target_dir}/${sub}/preprocessing_dwi/dwi_mt_off.bval";
    bvec_mt_off="${target_dir}/${sub}/preprocessing_dwi/dwi_mt_off.bvec";
    dwi_mt_on="${target_dir}/${sub}/preprocessing_dwi/dwi_mt_on_upsample.nii.gz";
    bval_mt_on="${target_dir}/${sub}/preprocessing_dwi/dwi_mt_off.bval";
    bvec_mt_on="${target_dir}/${sub}/preprocessing_dwi/dwi_mt_on.bvec";

    # ---------------------Bet-----------------------
    if [ ! -f "b0_mt_off_brain.nii.gz" ]; then
        echo "Extract b0";
        scil_dwi_extract_b0 $dwi_mt_off $bval_mt_off $bvec_mt_off b0_mt_off.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;

        echo "Bet b0";
        bet b0_mt_off.nii.gz b0_mt_off_brain -m;
    fi
    b0_mt_off="${target_dir}/${sub}/preprocessing_dwi/b0_mt_off_brain.nii.gz";
    brain_mask_mt_off="${target_dir}/${sub}/preprocessing_dwi/b0_mt_off_brain_mask.nii.gz";

    # --------------------Re-stride and cleanup-----------------------
    cd ${target_dir}/${sub};
    mkdir -p dwi;
    cd ${target_dir}/${sub}/dwi;
    if [ ! -f "dwi_mt_off.nii.gz" ]; then
        echo "Re-stride";
        # MT-off
        dwi_off=$dwi_mt_off;
        bval_off=$bval_mt_off;
        bvec_off=$bvec_mt_off;
        mask_off=$brain_mask_mt_off;
        b0=$b0_mt_off;
        mrconvert -strides 1,2,3,4 $dwi_off dwi_mt_off.nii.gz;
        mrconvert -strides 1,2,3 $mask_off b0_brain_mask.nii.gz;
        mrconvert -strides 1,2,3 $b0 b0.nii.gz;
        scil_gradients_modify_axes $bvec_off dwi_mt_off.bvec -1 2 3; # in case of -1 2 3 4 strides
        cp $bval_off dwi_mt_off.bval;
        # Some DWI in the dataset have NaN values that need to be removed
		python $code_dir/python_scripts/remove_nans_from_dwi.py dwi_mt_off.nii.gz dwi_mt_off.nii.gz -f;
        # MT-on
        dwi_on=$dwi_mt_on;
        bval_on=$bval_mt_on;
        bvec_on=$bvec_mt_on;
        mrconvert -strides 1,2,3,4 $dwi_on dwi_mt_on.nii.gz;
        scil_gradients_modify_axes $bvec_on dwi_mt_on.bvec -1 2 3; # in case of -1 2 3 4 strides
        cp $bval_on dwi_mt_on.bval;
        # Some DWI in the dataset have NaN values that need to be removed
		python $code_dir/python_scripts/remove_nans_from_dwi.py dwi_mt_on.nii.gz dwi_mt_on.nii.gz -f;
    fi

done;
