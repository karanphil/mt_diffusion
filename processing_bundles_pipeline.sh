#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh
# and processing_rbx_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"
# The third argument is the rbx_flow directory to create and use (full path)
rbx_data_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/rbx_flow"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------Bundles PROCESSING-------------------------
echo "Bundles PROCESSING";
for sub in $subs; 
    do echo $sub;

    fa="${target_dir}/${sub}/dti/fa.nii.gz";
    tractogram="${target_dir}/${sub}/tractography/local_tracking.trk";
    peaks_mt_off="${target_dir}/${sub}/fodf_metrics_mt_off/peaks.nii.gz";

    # ----------------------Get bundles------------------------
    cd ${target_dir}/${sub};
    mkdir -p bundles;
    cd ${target_dir}/${sub}/bundles;
    if [ ! -f "AF_L.trk" ]; then
        echo "Copy recognize bundles from rbx_flow output";
        cp -L ${rbx_data_dir}/output/results_rbx/${sub}/Recognize_Bundles/*.trk .;
    fi

    # -----------------Add dps to bundles and go to 2mm iso-------------------
    if [ ! -d "removed_bundles" ]; then
        echo "Add dps to bundles and resample to 2mm iso";
        bundles=$(ls *.trk);
        python ${code_dir}/python_scripts/add_dps_to_bundle.py $tractogram ${rbx_data_dir}/output/results_rbx/${sub}/Recognize_Bundles/results.json . --in_bundles $bundles -f;
        for b in $bundles;
            do bundle_name=${b%".trk"};
            # Compute outliers rejection
            scil_bundle_reject_outliers.py $b $b --alpha 0.5 -f;
            # Export dps files
            scil_tractogram_dps_math.py $b export sift2 --out_dps_file ${bundle_name}_sift2_weights.txt -f;
            # Resave bundles with reference in 2mm iso
            n=${bundle_name}.tck;
            scil_tractogram_convert.py $b $n;
            scil_tractogram_convert.py $n $b --reference $fa -f;
            rm *.tck;
            # Add SIFT2 weights
            scil_tractogram_dps_math.py $b import sift2 --out_tractogram $b --in_dps_file ${bundle_name}_sift2_weights.txt -f;
            rm ${bundle_name}_sift2_weights.txt;

        done;
        mkdir -p removed_bundles;
        # Remove CR bundle from the main folder, as it is not used in the analyses
        mv CR_* removed_bundles/;
    fi

    # ----------------------Fixel analysis------------------------
    cd ${target_dir}/${sub};
    mkdir -p fixel_analysis;
    if [ ! -f "fixel_analysis/single_bundle_mask_voxel-norm_WM.nii.gz" ]; then
        echo "Fixel analysis";
        # These thresholds ensure that fixels populated only by a few streamlines of a given bundle are not counted for the masks. These do not affect the maps.
        rel_thr=0.1;
        abs_thr=1.5;
        scil_bundle_fixel_analysis.py $peaks_mt_off --in_bundles ${target_dir}/${sub}/bundles/*.trk --processes 8 --single_bundle --split_bundles --rel_thr 0.1 --abs_thr 1.5 --norm voxel none --out_dir fixel_analysis/ -f --dps_key sift2;
        cd ${target_dir}/${sub}/fixel_analysis;
        cp single_bundle_mask_voxel-norm_WM.nii.gz tmp1.nii.gz;
        cp single_bundle_mask_none-norm_WM.nii.gz tmp2.nii.gz;
        rm single_bundle_*.nii.gz;
        mv tmp1.nii.gz single_bundle_mask_voxel-norm_WM.nii.gz;
        mv tmp2.nii.gz single_bundle_mask_none-norm_WM.nii.gz;
        rm fixel_density_map_none-norm_*.nii.gz;
        rm fixel_density_mask_none-norm_*.nii.gz;
        rm voxel_density_map_none-norm_*.nii.gz;
        rm voxel_density_mask_none-norm_*.nii.gz;
    fi

done;
