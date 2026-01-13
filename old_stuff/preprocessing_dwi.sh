#!/bin/bash
# À rouler dans /home/local/USHERBROOKE/karp2601/Samsung/data/mt-diff-mcgill

b0_thr_extract_b0=10;

processed_data_dir="processed_data/preprocessing_dwi";
cd $processed_data_dir;
renamed_data_dir="../../renamed_data";

dwi_mt_off="${renamed_data_dir}/mt_off_dwi.nii.gz";
dwi_mt_on="${renamed_data_dir}/mt_on_dwi.nii.gz";
bval="${renamed_data_dir}/mt_off_dwi.bval";
bvec="${renamed_data_dir}/mt_off_dwi.bvec";
rev_b0="${renamed_data_dir}/mt_off_revb0.nii.gz";

# Merge for denoising
if [ ! -f "mt_off_dwi_dn.nii.gz" ]; then
echo "Denoising";
mrcat -axis 3 $dwi_mt_off $dwi_mt_on dwi_merge.nii.gz -f;
dwidenoise dwi_merge.nii.gz dwi_merge_dn.nii.gz -f;
rm dwi_merge.nii.gz;
mrconvert -coord 3 0:30 dwi_merge_dn.nii.gz mt_off_dwi_dn.nii.gz -force;
mrconvert -coord 3 31:end dwi_merge_dn.nii.gz mt_on_dwi_dn.nii.gz -force;
rm dwi_merge_dn.nii.gz;
fi
dwi_mt_off="mt_off_dwi_dn.nii.gz";
dwi_mt_on="mt_on_dwi_dn.nii.gz";

# dwidenoise $rev_b0 mt_off_revb0_dn.nii.gz -f;
# rev_b0="mt_off_revb0_dn.nii.gz";

echo "Extract b0";
scil_dwi_extract_b0.py $dwi_mt_off $bval $bvec mt_off_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
b0_mt_off="mt_off_b0.nii.gz";

# Add b0 MT-off to MT-on DWI
scil_volume_math.py concatenate $b0_mt_off $dwi_mt_on mt_on_dwi_extended.nii.gz  -f;
sed -e 's/^/0 /' $bvec > mt_on_dwi_expended.bvec;
sed -e 's/^/0 /' $bval > mt_on_dwi_expended.bval;
dwi_mt_on_ext="mt_on_dwi_extended.nii.gz";
bvec_ext="mt_on_dwi_expended.bvec";
bval_ext="mt_on_dwi_expended.bval";

echo "Bet b0";
bet $b0_mt_off b0_mt_off_brain -m;
brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

# For the moment, the images seem very well aligned, so skip registration.
echo "Topup";
scil_dwi_prepare_topup_command.py $b0_mt_off $rev_b0 --out_script --out_prefix topup -f;
sh topup.sh;

# Run eddy on MT-off
scil_dwi_prepare_eddy_command.py $dwi_mt_off $bval $bvec $brain_mask_mt_off --eddy_cmd eddy_cpu\
            --b0_thr $b0_thr_extract_b0\
            --out_script --fix_seed\
            --lsr_resampling --slice_drop_correction\
            --topup topup\
            --out_prefix dwi_mt_off -f;
sh eddy.sh;

cp $bval dwi_mt_off.bval;
cp dwi_mt_off.eddy_rotated_bvecs dwi_mt_off.bvec;
dwi_mt_off="dwi_mt_off.nii.gz";

echo "Extract b0";
scil_dwi_extract_b0.py $dwi_mt_off dwi_mt_off.bval dwi_mt_off.bvec mt_off_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
b0_mt_off="mt_off_b0.nii.gz";

echo "Bet b0";
bet $b0_mt_off b0_mt_off_brain -m;
brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

# Resample dwi-mt-off
scil_volume_resample.py $dwi_mt_off dwi_mt_off_upsample.nii.gz --voxel_size 1 -f;

# # Normalizing DWI MT-off
# mrcalc dwi_mt_off_upsample.nii.gz b0_mt_off_brain.nii.gz -div dwi_mt_off_norm.nii.gz -force;
# mrcalc dwi_mt_off_norm.nii.gz b0_mt_off_brain_mask.nii.gz -mult dwi_mt_off_norm.nii.gz -force;
# scil_volume_math.py upper_clip dwi_mt_off_norm.nii.gz 1 dwi_mt_off_norm.nii.gz -f;
# scil_volume_math.py lower_clip dwi_mt_off_norm.nii.gz 0 dwi_mt_off_norm.nii.gz -f;

# Run eddy on MT-on
scil_dwi_prepare_eddy_command.py $dwi_mt_on_ext $bval_ext $bvec_ext $brain_mask_mt_off --eddy_cmd eddy_cpu\
            --b0_thr $b0_thr_extract_b0\
            --out_script --fix_seed\
            --lsr_resampling --slice_drop_correction\
            --topup topup\
            --out_prefix dwi_mt_on -f;
sh eddy.sh;

cp dwi_mt_off.bval dwi_mt_on.bval;
sed -E 's/^[[:space:]]*[-+]?[0-9]+(\.[0-9]*)?[eE][-+]?[0-9]+[[:space:]]*//' dwi_mt_on.eddy_rotated_bvecs > dwi_mt_on.bvec;
mrconvert -coord 3 1:31 dwi_mt_on.nii.gz dwi_mt_on.nii.gz -force;

# MT-OFF
echo "Extract b0";
scil_dwi_extract_b0.py dwi_mt_off_upsample.nii.gz dwi_mt_off.bval dwi_mt_off.bvec mt_off_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
b0_mt_off="mt_off_b0.nii.gz";
echo "Bet b0";
bet $b0_mt_off b0_mt_off_brain -m;
brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

# Resample dwi-mt-on
scil_volume_resample.py dwi_mt_on.nii.gz dwi_mt_on_upsample.nii.gz --voxel_size 1 -f;

# # Normalizing DWI MT-on
# mrcalc dwi_mt_on_upsample.nii.gz b0_mt_off_brain.nii.gz -div dwi_mt_on_norm.nii.gz -force;
# mrcalc dwi_mt_on_norm.nii.gz b0_mt_off_brain_mask.nii.gz -mult dwi_mt_on_norm.nii.gz -force;
# scil_volume_math.py upper_clip dwi_mt_on_norm.nii.gz 1 dwi_mt_on_norm.nii.gz -f;
# scil_volume_math.py lower_clip dwi_mt_on_norm.nii.gz 0 dwi_mt_on_norm.nii.gz -f;

# Resampling sphere of DWI MT-on to the MT-off bvecs to take Eddy into acount and allow for substraction.
scil_dwi_extract_b0.py dwi_mt_on_upsample.nii.gz dwi_mt_on.bval dwi_mt_on.bvec mt_on_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
scil_dwi_to_sh.py dwi_mt_on_upsample.nii.gz dwi_mt_on.bval dwi_mt_on.bvec dwi_mt_on_sh.nii.gz --sh_order 6 --smooth 0 --mask b0_mt_off_brain_mask.nii.gz -f;
scil_sh_to_sf.py dwi_mt_on_sh.nii.gz dwi_mt_on_upsample_resample.nii.gz --in_bvec dwi_mt_off.bvec --in_bval dwi_mt_off.bval --out_bval toto.bval --in_b0 mt_on_b0.nii.gz --processes 8 -f;
rm toto.bval;

scil_volume_math.py subtraction dwi_mt_off_upsample.nii.gz dwi_mt_on_upsample_resample.nii.gz dwi_diff.nii.gz -f;
mrcalc dwi_diff.nii.gz dwi_mt_off_upsample.nii.gz -div dwi_MTR.nii.gz -force;
mrcalc dwi_MTR.nii.gz b0_mt_off_brain_mask.nii.gz -mult dwi_MTR.nii.gz -force;
scil_volume_math.py lower_clip dwi_MTR.nii.gz 0 dwi_MTR.nii.gz -f;
scil_volume_math.py upper_clip dwi_MTR.nii.gz 1 dwi_MTR.nii.gz -f;

scil_dwi_extract_b0.py dwi_MTR.nii.gz dwi_mt_off.bval dwi_mt_off.bvec MTR_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
# scil_volume_math.py upper_threshold_eq MTR_b0.nii.gz 0.2 MTR_mask.nii.gz -f;
# scil_volume_math.py invert MTR_mask.nii.gz MTR_mask.nii.gz -f;
# scil_volume_math.py multiplication $brain_mask_mt_off MTR_mask.nii.gz wm_mask.nii.gz --data_type uint8 -f;