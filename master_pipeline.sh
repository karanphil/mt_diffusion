#!/bin/bash
# À rouler dans /home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing

target_dir="/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing";
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------DIFFUSION PREPROCESSING--------------------------
b0_thr_extract_b0=10;

echo "DIFFUSION PREPROCESSING";
for sub in $subs; 
    do echo $sub;
    # cd ${target_dir}/${sub};
    # mkdir preprocessing_dwi;
    # cd ${target_dir}/${sub}/preprocessing_dwi;
    # renamed_data_dir="../renamed_data";

    # dwi_mt_off="${renamed_data_dir}/mt_off_dwi.nii.gz";
    # dwi_mt_on="${renamed_data_dir}/mt_on_dwi.nii.gz";
    # bval="${renamed_data_dir}/mt_off_dwi.bval";
    # bvec="${renamed_data_dir}/mt_off_dwi.bvec";
    # rev_b0="${renamed_data_dir}/mt_off_revb0.nii.gz";

    # # Merge for denoising
    # if [ ! -f "mt_off_dwi_dn.nii.gz" ]; then
    # echo "Denoising";
    # mrcat -axis 3 $dwi_mt_off $dwi_mt_on dwi_merge.nii.gz -f;
    # dwidenoise dwi_merge.nii.gz dwi_merge_dn.nii.gz -f;
    # rm dwi_merge.nii.gz;
    # mrconvert -coord 3 0:30 dwi_merge_dn.nii.gz mt_off_dwi_dn.nii.gz -force;
    # mrconvert -coord 3 31:end dwi_merge_dn.nii.gz mt_on_dwi_dn.nii.gz -force;
    # rm dwi_merge_dn.nii.gz;
    # fi
    # dwi_mt_off="mt_off_dwi_dn.nii.gz";
    # dwi_mt_on="mt_on_dwi_dn.nii.gz";

    # echo "Extract b0";
    # scil_dwi_extract_b0.py $dwi_mt_off $bval $bvec mt_off_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
    # b0_mt_off="mt_off_b0.nii.gz";

    # # Add b0 MT-off to MT-on DWI
    # scil_volume_math.py concatenate $b0_mt_off $dwi_mt_on mt_on_dwi_extended.nii.gz  -f;
    # sed -e 's/^/0 /' $bvec > mt_on_dwi_expended.bvec;
    # sed -e 's/^/0 /' $bval > mt_on_dwi_expended.bval;
    # dwi_mt_on_ext="mt_on_dwi_extended.nii.gz";
    # bvec_ext="mt_on_dwi_expended.bvec";
    # bval_ext="mt_on_dwi_expended.bval";

    # echo "Bet b0";
    # bet $b0_mt_off b0_mt_off_brain -m;
    # brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

    # # For the moment, the images seem very well aligned, so skip registration.
    # echo "Topup";
    # scil_dwi_prepare_topup_command.py $b0_mt_off $rev_b0 --out_script --out_prefix topup -f;
    # sh topup.sh;

    # # Run eddy on MT-off
    # echo "Eddy";
    # scil_dwi_prepare_eddy_command.py $dwi_mt_off $bval $bvec $brain_mask_mt_off --eddy_cmd eddy_cpu\
    #             --b0_thr $b0_thr_extract_b0\
    #             --out_script --fix_seed\
    #             --lsr_resampling --slice_drop_correction\
    #             --topup topup\
    #             --out_prefix dwi_mt_off -f;
    # sh eddy.sh;

    # cp $bval dwi_mt_off.bval;
    # cp dwi_mt_off.eddy_rotated_bvecs dwi_mt_off.bvec;
    # dwi_mt_off="dwi_mt_off.nii.gz";

    # echo "Extract b0";
    # scil_dwi_extract_b0.py $dwi_mt_off dwi_mt_off.bval dwi_mt_off.bvec mt_off_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
    # b0_mt_off="mt_off_b0.nii.gz";

    # echo "Bet b0";
    # bet $b0_mt_off b0_mt_off_brain -m;
    # brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

    # # Run eddy on MT-on
    # echo "Eddy";
    # scil_dwi_prepare_eddy_command.py $dwi_mt_on_ext $bval_ext $bvec_ext $brain_mask_mt_off --eddy_cmd eddy_cpu\
    #             --b0_thr $b0_thr_extract_b0\
    #             --out_script --fix_seed\
    #             --lsr_resampling --slice_drop_correction\
    #             --topup topup\
    #             --out_prefix dwi_mt_on -f;
    # sh eddy.sh;

    # cp dwi_mt_off.bval dwi_mt_on.bval;
    # sed -E 's/^[[:space:]]*[-+]?[0-9]+(\.[0-9]*)?[eE][-+]?[0-9]+[[:space:]]*//' dwi_mt_on.eddy_rotated_bvecs > dwi_mt_on.bvec;
    # mrconvert -coord 3 1:31 dwi_mt_on.nii.gz dwi_mt_on.nii.gz -force;

    # # Resampling sphere of DWI MT-on to the MT-off bvecs to take Eddy into acount and allow for substraction.
    # scil_dwi_extract_b0.py dwi_mt_on.nii.gz dwi_mt_on.bval dwi_mt_on.bvec mt_on_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
    # scil_dwi_to_sh.py dwi_mt_on.nii.gz dwi_mt_on.bval dwi_mt_on.bvec dwi_mt_on_sh.nii.gz --sh_order 6 --smooth 0 --mask b0_mt_off_brain_mask.nii.gz -f;
    # scil_sh_to_sf.py dwi_mt_on_sh.nii.gz dwi_mt_on_resample.nii.gz --in_bvec dwi_mt_off.bvec --in_bval dwi_mt_off.bval --out_bval toto.bval --in_b0 mt_on_b0.nii.gz --processes 8 -f;
    # rm toto.bval;
    # rm dwi_mt_on_sh.nii.gz;
    # rm mt_on_b0.nii.gz;

    # echo "Resample";
    # # Resample dwi-mt-off
    # scil_volume_resample.py $dwi_mt_off dwi_mt_off_upsample.nii.gz --voxel_size 2 -f;
    # # Resample dwi-mt-on
    # scil_volume_resample.py dwi_mt_on_resample.nii.gz dwi_mt_on_upsample.nii.gz --voxel_size 2 -f;

    # # MT-OFF
    # echo "Extract b0";
    # scil_dwi_extract_b0.py dwi_mt_off_upsample.nii.gz dwi_mt_off.bval dwi_mt_off.bvec mt_off_b0.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
    # b0_mt_off="mt_off_b0.nii.gz";
    # echo "Bet b0";
    # bet $b0_mt_off b0_mt_off_brain -m;
    # brain_mask_mt_off="b0_mt_off_brain_mask.nii.gz";

    # # Normalizing DWI MT-off
    # mrcalc dwi_mt_off_upsample.nii.gz b0_mt_off_brain.nii.gz -div dwi_mt_off_norm.nii.gz -force;
    # mrcalc dwi_mt_off_norm.nii.gz b0_mt_off_brain_mask.nii.gz -mult dwi_mt_off_norm.nii.gz -force;
    # scil_volume_math.py upper_clip dwi_mt_off_norm.nii.gz 1 dwi_mt_off_norm.nii.gz -f;
    # scil_volume_math.py lower_clip dwi_mt_off_norm.nii.gz 0 dwi_mt_off_norm.nii.gz -f;

    # # Normalizing DWI MT-on
    # mrcalc dwi_mt_on_upsample.nii.gz b0_mt_off_brain.nii.gz -div dwi_mt_on_norm.nii.gz -force;
    # mrcalc dwi_mt_on_norm.nii.gz b0_mt_off_brain_mask.nii.gz -mult dwi_mt_on_norm.nii.gz -force;
    # scil_volume_math.py upper_clip dwi_mt_on_norm.nii.gz 1 dwi_mt_on_norm.nii.gz -f;
    # scil_volume_math.py lower_clip dwi_mt_on_norm.nii.gz 0 dwi_mt_on_norm.nii.gz -f;

    # # Re-striding data
    # cd ${target_dir}/${sub};
    # mkdir dwi;
    # cd ${target_dir}/${sub}/dwi;
    # MT-off
    # dwi_off="../preprocessing_dwi/dwi_mt_off_norm.nii.gz";
    # bval_off="../preprocessing_dwi/dwi_mt_off.bval";
    # bvec_off="../preprocessing_dwi/dwi_mt_off.bvec";
    # mask_off="../preprocessing_dwi/b0_mt_off_brain_mask.nii.gz";
    # b0="../preprocessing_dwi/b0_mt_off_brain.nii.gz";
    # mrconvert -strides 1,2,3,4 $dwi_off dwi_mt_off.nii.gz;
    # mrconvert -strides 1,2,3 $mask_off b0_brain_mask.nii.gz;
    # mrconvert -strides 1,2,3 $b0 b0.nii.gz;
    # scil_gradients_modify_axes.py $bvec_off dwi.bvec -1 2 3; # in case of -1 2 3 4 strides
    # cp $bval_off dwi.bval;
    # # MT-on
    # dwi_on="../preprocessing_dwi/dwi_mt_on_norm.nii.gz";
    # mrconvert -strides 1,2,3,4 $dwi_on dwi_mt_on.nii.gz;

    dwi_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.nii.gz";
    dwi_mt_on="${target_dir}/${sub}/dwi/dwi_mt_on.nii.gz";
    bval="${target_dir}/${sub}/dwi/dwi.bval";
    bvec="${target_dir}/${sub}/dwi/dwi.bvec";
    mask="${target_dir}/${sub}/dwi/b0_brain_mask.nii.gz";
    b0="${target_dir}/${sub}/dwi/b0.nii.gz";

    # # Computing PA stuff
    # cd ${target_dir}/${sub};
    # mkdir powder_average;
    # cd ${target_dir}/${sub}/powder_average;
    # scil_volume_math.py subtraction $dwi_mt_off $dwi_mt_on powder_averaged_mtr.nii.gz --data_type float32 -f;
    # mrcalc powder_averaged_mtr.nii.gz $dwi_mt_off -div powder_averaged_mtr.nii.gz -force;
    # mrcalc powder_averaged_mtr.nii.gz $mask -mult powder_averaged_mtr.nii.gz -force;
    # scil_dwi_extract_b0.py powder_averaged_mtr.nii.gz $bval $bvec b0_mtr.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
    # scil_volume_math.py lower_clip b0_mtr.nii.gz 0 b0_mtr.nii.gz -f;
    # scil_volume_math.py upper_clip b0_mtr.nii.gz 1 b0_mtr.nii.gz -f;
    # scil_dwi_powder_average.py powder_averaged_mtr.nii.gz $bval powder_averaged_mtr.nii.gz --mask $mask -f;
    # scil_volume_math.py lower_clip powder_averaged_mtr.nii.gz 0 powder_averaged_mtr.nii.gz -f;
    # scil_volume_math.py upper_clip powder_averaged_mtr.nii.gz 1 powder_averaged_mtr.nii.gz -f;

    # # Resample DWI mt-off to 1mm iso for tractography
    # echo "Resample";
    # cd ${target_dir}/${sub};
    # mkdir dwi_for_tractography;
    # cd ${target_dir}/${sub}/dwi_for_tractography;
    # # Resample dwi-mt-off
    # mrconvert -strides 1,2,3,4 ../preprocessing_dwi/dwi_mt_off.nii.gz dwi_mt_off.nii.gz;
    # scil_volume_resample.py dwi_mt_off.nii.gz dwi_mt_off_upsample_for_tractography.nii.gz --voxel_size 1 -f;
    # echo "Extract b0";
    # scil_dwi_extract_b0.py dwi_mt_off_upsample_for_tractography.nii.gz $bval $bvec mt_off_b0_for_tractography.nii.gz --mean --b0_threshold $b0_thr_extract_b0 --skip_b0_check;
    # b0_mt_off="mt_off_b0_for_tractography.nii.gz";
    # echo "Bet b0";
    # bet $b0_mt_off b0_mt_off_for_tractography_brain -m;
    # brain_mask_mt_off="b0_mt_off_for_tractography_brain_mask.nii.gz";
    dwi_for_tractography="${target_dir}/${sub}/dwi_for_tractography/dwi_mt_off_upsample_for_tractography.nii.gz";
    mask_for_tractography="${target_dir}/${sub}/dwi_for_tractography/b0_mt_off_for_tractography_brain_mask.nii.gz";
    b0_for_tractography="${target_dir}/${sub}/dwi_for_tractography/b0_mt_off_for_tractography_brain.nii.gz";

    # !!! Some DTI will crash because of NaNs inside the mask. !!!
    # This is then needed in python:
    # '
    # import nibabel as nib
    # import numpy as np
    # vol = nib.load("dwi_mt_off_norm.nii.gz")
    # dwi = vol.get_fdata()
    # mask_vol = nib.load("b0_mt_off_brain_mask.nii.gz")
    # mask = mask_vol.get_fdata().astype(bool)
    # new_dwi = np.where(np.isnan(dwi), 0, dwi)
    # nib.save(nib.Nifti1Image(new_dwi, vol.affine), "dwi_mt_off_norm.nii.gz")
    # '

    # # Compute DTI
    # cd ${target_dir}/${sub};
    # mkdir dti;
    # cd ${target_dir}/${sub}/dti;
    # scil_dti_metrics.py $dwi_mt_off $bval $bvec --mask $mask -f --not_all --fa fa.nii.gz --md md.nii.gz --rgb rgb.nii.gz;
    fa="${target_dir}/${sub}/dti/fa.nii.gz";
    md="${target_dir}/${sub}/dti/md.nii.gz";

    # # Compute DTI for tractography
    # cd ${target_dir}/${sub};
    # mkdir dti_for_tractography;
    # cd ${target_dir}/${sub}/dti_for_tractography;
    # scil_dti_metrics.py $dwi_for_tracto $bval $bvec --mask $mask_for_tracto -f --not_all --fa fa.nii.gz --md md.nii.gz --rgb rgb.nii.gz;
    fa_for_tractography="${target_dir}/${sub}/dti_for_tractography/fa.nii.gz";
    md_for_tractography="${target_dir}/${sub}/dti_for_tractography/md.nii.gz";

    # echo "T1";
    t1="${target_dir}/${sub}/renamed_data/t1.nii.gz";

    # cd ${target_dir}/${sub};
    # mkdir preprocessing_t1;
    # cd ${target_dir}/${sub}/preprocessing_t1;

    # cd ../../..;
    # singularity exec -B /data/karp2601/stockage/mt-diff-mcgill/ ~/Research/containers/scilus_2.0.2_from_docker.sif bash code/mt_diffusion/t1_pipeline.sh $t1 $fa $b0 $mask $sub preprocessing_t1;
    t1="${target_dir}/${sub}/preprocessing_t1/register_natif/outputWarped.nii.gz";
    wm_mask="${target_dir}/${sub}/preprocessing_t1/register_natif/wm_mask.nii.gz";

    # cd ${target_dir}/${sub};
    # mkdir preprocessing_t1_for_tractography;
    # cd ${target_dir}/${sub}/preprocessing_t1_for_tractography;

    # cd ../../..;
    # singularity exec -B /data/karp2601/stockage/mt-diff-mcgill/ ~/Research/containers/scilus_2.0.2_from_docker.sif bash code/mt_diffusion/t1_pipeline.sh $t1 $fa_for_tractography $b0_for_tractography $mask_for_tractography $sub preprocessing_t1_for_tractography;
    wm_mask_for_tractography="${target_dir}/${sub}/preprocessing_t1_for_tractography/register_natif/wm_mask.nii.gz";

    # # Compute CSD
    # echo "CSD";
    # cd ${target_dir}/${sub};
    # mkdir fodf;
    # cd ${target_dir}/${sub}/fodf;
    # scil_frf_ssst.py $dwi_mt_off $bval $bvec frf.txt --mask $mask --mask_wm $wm_mask --roi_radii 15 15 10 -f;
    # scil_fodf_ssst.py $dwi_mt_off $bval $bvec frf.txt fodf_mt_off.nii.gz --mask $mask --processes 8 --sh_order 6 -f;
    # scil_fodf_ssst.py $dwi_mt_on $bval $bvec frf.txt fodf_mt_on.nii.gz --mask $mask --processes 8 --sh_order 6 -f;
    fodf_mt_off="${target_dir}/${sub}/fodf/fodf_mt_off.nii.gz";
    fodf_mt_on="${target_dir}/${sub}/fodf/fodf_mt_on.nii.gz";

    # # Compute FODF metrics for fixel analysis
    # echo "FODF metrics for fixel analysis";
    # cd ${target_dir}/${sub};
    # mkdir fodf_metrics_mt_off;
    # cd ${target_dir}/${sub}/fodf_metrics_mt_off;
    # scil_fodf_max_in_ventricles $fodf_mt_off $fa $md --md_threshold 0.0025 --max_value_output max_fodf_in_ventricles.txt --in_mask ${target_dir}/${sub}/preprocessing_t1/register_natif/csf_mask.nii.gz --use_median -f;
    # max_value=$(cat max_fodf_in_ventricles.txt);    
    # a_threshold=$(echo 2*${max_value}|bc);
    # scil_fodf_metrics $fodf_mt_off --mask $mask --abs_peaks_and_values --at $a_threshold -f --processes 8;
    peaks_mt_off="${target_dir}/${sub}/fodf_metrics_mt_off/peaks.nii.gz";

    # # Compute FODF metrics of mt_on
    # echo "FODF metrics of mt_on";
    # cd ${target_dir}/${sub};
    # mkdir fodf_metrics_mt_on;
    # cd ${target_dir}/${sub}/fodf_metrics_mt_on;
    # scil_fodf_max_in_ventricles $fodf_mt_on $fa $md --md_threshold 0.0025 --max_value_output max_fodf_in_ventricles.txt --in_mask ${target_dir}/${sub}/preprocessing_t1/register_natif/csf_mask.nii.gz --use_median -f;
    # max_value=$(cat max_fodf_in_ventricles.txt);    
    # a_threshold=$(echo 2*${max_value}|bc);
    # # !!! Needed to add a check for nans in scil_fodf_metrics.py because some images have nans in the fodf.
    # scil_fodf_metrics $fodf_mt_on --mask $mask --abs_peaks_and_values --at $a_threshold -f --processes 8;
    peaks_mt_on="${target_dir}/${sub}/fodf_metrics_mt_on/peaks.nii.gz";

    # # Compute CSD for tractography
    # echo "CSD for tractography";
    # cd ${target_dir}/${sub};
    # mkdir fodf_for_tractography;
    # cd ${target_dir}/${sub}/fodf_for_tractography;
    # scil_frf_ssst.py $dwi_for_tractography $bval $bvec frf.txt --mask $mask_for_tractography --mask_wm $wm_mask_for_tractography --roi_radii 30 30 20 -f;
    # scil_fodf_ssst.py $dwi_for_tractography $bval $bvec frf.txt fodf.nii.gz --mask $mask_for_tractography --processes 8 --sh_order 6 -f;
    fodf_for_tractography="${target_dir}/${sub}/fodf_for_tractography/fodf.nii.gz";

    # # Compute FODF metrics for tractography
    # echo "FODF metrics for tractography";
    # cd ${target_dir}/${sub};
    # mkdir fodf_metrics_for_tractography;
    # cd ${target_dir}/${sub}/fodf_metrics_for_tractography;
    # scil_fodf_max_in_ventricles $fodf_for_tractography $fa_for_tractography $md_for_tractography --md_threshold 0.0025 --max_value_output max_fodf_in_ventricles.txt --in_mask ${target_dir}/${sub}/preprocessing_t1_for_tractography/register_natif/csf_mask.nii.gz --use_median -f;
    # max_value=$(cat max_fodf_in_ventricles.txt);    
    # a_threshold=$(echo 2*${max_value}|bc);
    # scil_fodf_metrics $fodf_for_tractography --mask $mask_for_tractography --abs_peaks_and_values --at $a_threshold -f --processes 8;

    # # Tractography
    # echo "Tractography";
    # cd ${target_dir}/${sub};
    # mkdir tractography;
    # cd ${target_dir}/${sub}/tractography;
    # # scil_tracking_pft_maps $wm_mask_for_tractography ${target_dir}/${sub}/preprocessing_t1_for_tractography/register_natif/gm_mask.nii.gz ${target_dir}/${sub}/preprocessing_t1_for_tractography/register_natif/csf_mask.nii.gz --include map_include.nii.gz --exclude map_exclude.nii.gz --interface interface.nii.gz -f;
    # # scil_volume_math convert $wm_mask_for_tractography pft_seeding_mask.nii.gz --data_type uint8 -f;
    # # scil_volume_math union pft_seeding_mask.nii.gz interface.nii.gz pft_seeding_mask.nii.gz --data_type uint8 -f;
    # # scil_tracking_pft $fodf_for_tractography pft_seeding_mask.nii.gz map_include.nii.gz map_exclude.nii.gz pft_tracking.trk --algo prob --npv 10 --seed 0 --step 0.5 --theta 20 --min_length 20 --max_length 200 --particles 15 --back 2 --forward 1 --compress 0.2 --sh_basis descoteaux07 -f;
    # # scil_tractogram_remove_invalid pft_tracking.trk pft_tracking.trk --remove_single_point -f;
    # scil_tracking_local $fodf_for_tractography $wm_mask_for_tractography $wm_mask_for_tractography local_tracking.trk --use_gpu --npv 10 -f;
    # scil_tractogram_remove_invalid local_tracking.trk local_tracking.trk --remove_single_point -f;
    tractogram="${target_dir}/${sub}/tractography/local_tracking.trk";

    # # Compute SIFT2
    # echo "SIFT2";
    # cd ${target_dir}/${sub}/tractography;
    # fodf_tournier="${target_dir}/${sub}/fodf_tournier.nii.gz";
    # scil_sh_convert $fodf_for_tractography $fodf_tournier descoteaux07_legacy tournier07 -f;
    # tractogram_tck="${target_dir}/${sub}/local_tracking.tck";
    # scil_tractogram_convert $tractogram $tractogram_tck -f;
    # tcksift2 $tractogram_tck $fodf_tournier sift2_weights.txt -force;
    # scil_tractogram_dps_math $tractogram import "sift2" --in_dps_file sift2_weights.txt --out_tractogram $tractogram -f;
    # rm $tractogram_tck $fodf_tournier;

    # # Run rbx_flow on the side with nextflow.
    # # In ~/data/stockage/mt-diff-mcgill/rbx_flow
    # # nextflow run ~/Research/source/rbx_flow/main.nf --input ../input --atlas_directory ~/data/stockage/atlas_old -with-singularity /home/local/USHERBROOKE/karp2601/Research/containers/scilus_2.1.0.sif --register_processes 8  --rbx_processes  8
    # # Copy output to bundles folder.
    # echo "Copy recognize bundles from rbx_flow output";
    # cd ${target_dir}/${sub};
    # mkdir bundles;
    # cd ${target_dir}/${sub}/bundles;
    # rm -r *;
    # cp -L ~/data/stockage/mt-diff-mcgill/rbx_flow/output/results_rbx/${sub}/Recognize_Bundles/*.trk .;
    # bundles=$(ls);
    # python ../../../code/mt_diffusion/add_dps_to_bundle.py $tractogram ~/data/stockage/mt-diff-mcgill/rbx_flow/output/results_rbx/${sub}/Recognize_Bundles/results.json . --in_bundles $bundles -v -f;
    # for b in $bundles;
    #     do bundle_name=${b%".trk"};
    #     echo $bundle_name;
    #     # echo "Copy cleaned bundles";
    #     # rm -r $b;
    #     # cp -L ~/data/stockage/mt-diff-mcgill/rbx_flow/output/results_rbx/${sub}/Clean_Bundles/${sub}__${bundle_name}_cleaned.trk ${bundle_name}.trk;
    #     # scil_tractogram_math intersection $b ~/data/stockage/mt-diff-mcgill/rbx_flow/output/results_rbx/${sub}/Clean_Bundles/${sub}__${bundle_name}_cleaned.trk $b -p 3 -f -v;
    #     echo "Compute outliers rejection";
    #     scil_bundle_reject_outliers $b $b --alpha 0.5 -f;
    #     echo "Export dps files";
    #     scil_tractogram_dps_math $b export sift2 --out_dps_file ${bundle_name}_sift2_weights.txt -f;
    #     echo "Resave bundles with reference";
    #     n=${bundle_name}.tck;
    #     scil_tractogram_convert $b $n;
    #     scil_tractogram_convert $n $b --reference $fa -f;
    #     rm *.tck;
    #     echo "Add SIFT2 weights";
    #     scil_tractogram_dps_math $b import sift2 --out_tractogram $b --in_dps_file ${bundle_name}_sift2_weights.txt -f;
    #     rm ${bundle_name}_sift2_weights.txt;

    # done;
    # mkdir removed_bundles;
    # mv CR_* removed_bundles/;

    # # Fixel analysis
    # echo "Fixel analysis";
    # cd ${target_dir}/${sub};
    # mkdir fixel_analysis;
    # # cd ${target_dir}/${sub}/fixel_analysis;
    # scil_bundle_fixel_analysis $peaks_mt_off --in_bundles ${target_dir}/${sub}/bundles/*.trk  --processes 8 --single_bundle --split_bundles --rel_thr 0.1 --abs_thr 1.5 --norm voxel none --out_dir fixel_analysis/ -f --dps_key sift2;
    # cd ${target_dir}/${sub}/fixel_analysis;
    # cp single_bundle_mask_voxel-norm_WM.nii.gz tmp1.nii.gz;
    # cp single_bundle_mask_none-norm_WM.nii.gz tmp2.nii.gz;
    # rm single_bundle_*.nii.gz;
    # mv tmp1.nii.gz single_bundle_mask_voxel-norm_WM.nii.gz;
    # mv tmp2.nii.gz single_bundle_mask_none-norm_WM.nii.gz;
    # rm fixel_density_map_none-norm_*.nii.gz;
    # rm fixel_density_mask_none-norm_*.nii.gz;
    # rm voxel_density_map_none-norm_*.nii.gz;
    # rm voxel_density_map_voxel-norm_*.nii.gz;
    # rm voxel_density_mask_none-norm_*.nii.gz;
    # rm voxel_density_mask_voxel-norm_*.nii.gz;

    # Compute FODF MTR
    echo "Compute FODF MTR";
    cd ${target_dir}/${sub};
    mkdir mtr;
    cd ${target_dir}/${sub}/mtr;
    # rel_thr=0.1;
    # abs_thr=1.5;
    rel_thr=0.01;
    abs_thr=0.0;
    python ../../../code/mt_diffusion/compute_odf_mtr.py $fodf_mt_off $fodf_mt_on $peaks_mt_off $peaks_mt_on ${target_dir}/${sub}/fixel_analysis/fixel_density_maps_voxel-norm.nii.gz ${target_dir}/${sub}/fixel_analysis/fixel_density_maps_none-norm.nii.gz mtr_fodf.nii.gz mtr_peak_values.nii.gz mtr_peaks.nii.gz --mask $mask --rel_thr $rel_thr --abs_thr $abs_thr --min_angle 10 -f;

    # Compute bundle MTR
    echo "Compute bundle MTR";
    cd ${target_dir}/${sub}/bundles;
    bundles=$(ls *.trk);
    cd ${target_dir}/${sub}/mtr;
    for b in $bundles;
        do bundle_name=${b%".trk"};
        echo $bundle_name;
        python ../../../code/mt_diffusion/compute_bundle_mtr.py mtr_peak_values.nii.gz ${target_dir}/${sub}/fixel_analysis/fixel_density_mask_voxel-norm_${bundle_name}.nii.gz mtr_${bundle_name}.nii.gz;

    done;

    # # !!!!!!!!!!!!!!!!!! A rerouler avec nufo.nii.gz from fodf_metrics_mtr!!!!!!!!!!!!!!!!!!
    # # Clean crossing mask
    # echo "Clean crossing mask";
    # cd ${target_dir}/${sub};
    # mkdir clean_crossing_mask;
    # cd ${target_dir}/${sub}/clean_crossing_mask;
    # python ../../../code/mt_diffusion/compute_clean_crossing_mask.py ${target_dir}/${sub}/fodf_metrics_mt_off/nufo.nii.gz ${target_dir}/${sub}/fixel_analysis/nb_bundles_per_voxel_voxel-norm.nii.gz ${target_dir}/${sub}/fixel_analysis/fixel_density_maps_voxel-norm.nii.gz clean_crossing_mask.nii.gz --thr 0.7;

    # cd ../..;

done;
