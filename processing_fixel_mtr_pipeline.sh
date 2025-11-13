#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh,
# processing_rbx_pipeline.sh and processing_bundles_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------Fixel-MTR PROCESSING-------------------------
echo "Bundles PROCESSING";
for sub in $subs; 
    do echo $sub;

    dwi_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.nii.gz";
    bval_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.bval";
    bvec_mt_off="${target_dir}/${sub}/dwi/dwi_mt_off.bvec";
    fodf_mt_off="${target_dir}/${sub}/fodf/fodf_mt_off.nii.gz";
    fodf_mt_on="${target_dir}/${sub}/fodf/fodf_mt_on.nii.gz";
    peaks_mt_off="${target_dir}/${sub}/fodf_metrics_mt_off/peaks.nii.gz";
    peaks_mt_on="${target_dir}/${sub}/fodf_metrics_mt_on/peaks.nii.gz";
    nufo="${target_dir}/${sub}/fodf_metrics_mt_off/nufo.nii.gz";
    mask="${target_dir}/${sub}/dwi/b0_brain_mask.nii.gz";

    cd ${target_dir}/${sub}/bundles;
    bundles=$(ls *.trk);

    # -----------------------Fixel-MTR------------------------
    cd ${target_dir}/${sub};
    mkdir -p fixel_mtr;
    cd ${target_dir}/${sub}/fixel_mtr;
    if [ ! -f "mtr_peak_values.nii.gz" ]; then
        echo "Compute fixel-MTR";
        # These thresholds let almost everything go through, so we have a lot of fixels to look at with visualization.
        # The filtering will occur in the next steps using the masks instead of the maps. The masks are already thresholded.
        rel_thr=0.01;
        abs_thr=0.0;
        python ${code_dir}/python_scripts/compute_fixel_mtr.py $fodf_mt_off $fodf_mt_on $peaks_mt_off $peaks_mt_on ${target_dir}/${sub}/fixel_analysis/fixel_density_maps_voxel-norm.nii.gz ${target_dir}/${sub}/fixel_analysis/fixel_density_maps_none-norm.nii.gz mtr_fodf.nii.gz mtr_peak_values.nii.gz mtr_peaks.nii.gz --mask $mask --rel_thr $rel_thr --abs_thr $abs_thr --min_angle 20 -f;
    fi
    mtr_peak_values="${target_dir}/${sub}/fixel_mtr/mtr_peak_values.nii.gz";

    # --------------------Bundle fixel-MTR----------------------
    if [ ! -f "fixel_mtr_AF_L.nii.gz" ]; then
        echo "Compute bundle fixel-MTR";
        for b in $bundles;
            do bundle_name=${b%".trk"};
            python ${code_dir}/python_scripts/compute_bundle_fixel_mtr.py $mtr_peak_values ${target_dir}/${sub}/fixel_analysis/fixel_density_mask_voxel-norm_${bundle_name}.nii.gz fixel_mtr_${bundle_name}.nii.gz -f;

        done;
    fi

    # ----------------------Powder-average MTR---------------------
    cd ${target_dir}/${sub};
    mkdir -p mtr;
    cd ${target_dir}/${sub}/mtr;
    if [ ! -f "powder_averaged_mtr.nii.gz" ]; then
        echo "Powder Average MTR Calculation";
        echo "scil_dwi_powder_average.py $dwi_mt_off $bval_mt_off dwi_mt_off_pa.nii.gz --mask $mask -f";
        scil_dwi_powder_average.py $dwi_mt_off $bval_mt_off dwi_mt_off_pa.nii.gz --mask $mask -f;
        scil_dwi_powder_average.py $dwi_mt_on $bval_mt_on dwi_mt_on_pa.nii.gz --mask $mask -f;
        scil_volume_math.py subtraction dwi_mt_off_pa.nii.gz dwi_mt_on_pa.nii.gz powder_averaged_mtr.nii.gz --data_type float32 -f;
        mrcalc powder_averaged_mtr.nii.gz dwi_mt_off_pa.nii.gz -div powder_averaged_mtr.nii.gz -force;
        mrcalc powder_averaged_mtr.nii.gz $mask -mult powder_averaged_mtr.nii.gz -force;
    fi
    pa_mtr="${target_dir}/${sub}/mtr/powder_averaged_mtr.nii.gz";

    # -----------------------Bundle MTR--------------------------
    if [ ! -f "mtr_AF_L.nii.gz" ]; then
        echo "Compute bundle MTR";
        python ${code_dir}/python_scripts/compute_bundle_mtr.py $pa_mtr ${target_dir}/${sub}/fixel_analysis/voxel_density_masks_voxel-norm.nii.gz ${target_dir}/${sub}/fixel_analysis/bundles_LUT.txt . -f;
    fi

    # ----------------------Compute fixel-MTR difference-----------------------
    cd ${target_dir}/${sub}/fixel_mtr;
    if [ ! -f "mtr_peak_diffs.nii.gz" ]; then
        echo "Compute fixel-MTR difference for two first peaks";
        python ${code_dir}/python_scripts/compare_fixel_mtr_in_voxel.py $mtr_peak_values mtr_peak_diffs.nii.gz  mtr_peak_diff_mask.nii.gz mtr_peak_no_diff_mask.nii.gz mtr_peak_crossing_mask.nii.gz --min_diff 0.02 --min_mtr 0.2 -f;
    fi

    # ----------------------AFD fixel-----------------------
    cd ${target_dir}/${sub};
    mkdir -p afd_fixel;
    cd ${target_dir}/${sub}/afd_fixel;
    if [ ! -f "afd_fixel_AF_L.nii.gz" ]; then
        echo "Compute AFD fixel per bundle";
        for b in $bundles;
            do bundle_name=${b%".trk"};
            scil_bundle_mean_fixel_afd.py ${target_dir}/${sub}/bundles/${b} $fodf_mt_off afd_fixel_${bundle_name}.nii.gz -f;

        done;
    fi

    # ---------------------Compute crossing mask-----------------------
    cd ${target_dir}/${sub};
    mkdir -p clean_crossing_mask;
    cd ${target_dir}/${sub}/clean_crossing_mask;
    if [ ! -f "afd_fixel_AF_L.nii.gz" ]; then
        echo "Compute clean crossing mask";
        python ${code_dir}/python_scripts/compute_clean_crossing_mask.py $nufo ${target_dir}/${sub}/fixel_analysis/nb_bundles_per_voxel_voxel-norm.nii.gz ${target_dir}/${sub}/fixel_analysis/fixel_density_maps_voxel-norm.nii.gz clean_crossing_mask.nii.gz --thr 0.7;
    fi

done;