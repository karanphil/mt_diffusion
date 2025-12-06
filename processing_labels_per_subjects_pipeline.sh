#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh,
# processing_rbx_pipeline.sh, processing_bundles_pipeline.sh,
# processing_fixel_mtr_pipeline.sh, processing_registration_pipeline.sh and
# processing_averages_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"
# The third argument is the register_flow directory (full path)
register_data_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/register_flow"
cd $target_dir;
subs=$(ls -d hc*);

nb_sections=20; # This can be changed if needed

cd ${register_data_dir}/output/results_registration/hc18/Trks_into_template_space;
bundles=$(ls *.trk);
for sub in $subs; 
    do echo $sub;
    # ----------------------------Labels PROCESSING-------------------------
    cd ${register_data_dir}/output;
    mkdir -p processing_registration/${sub}/labels;
    cd ${register_data_dir}/output/processing_registration/${sub}/labels;
    mkdir -p tmp;
    if [ ! -f "labels_AF_L.nii.gz" ]; then
        echo "Labels PROCESSING";
        for bundle in $bundles;
            do b=${bundle%"_to_template.trk"};
            echo $b;
            scil_bundle_label_map.py ${register_data_dir}/output/results_registration/${sub}/Trks_into_template_space/${bundle} ${target_dir}/average_all_scans/centroid_${b}.trk tmp --nb_pts $nb_sections -f;
            cp tmp/labels_map.nii.gz labels_${b}.nii.gz;
        
        done;
    fi
    rm -r tmp;

    # -------------------------Bundles masks PROCESSING------------------------
    cd ${register_data_dir}/output;
    mkdir -p processing_registration/${sub}/bundles_masks;
    cd ${register_data_dir}/output/processing_registration/${sub}/bundles_masks;
    erosions="3 3 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 2 2 2 2 2 1 1 2 2 2 2 1 1" # This was hardcoded according to the size of the bundles.
    if [ ! -f "AF_L_mask_eroded.nii.gz" ]; then
        echo "Bundles masks PROCESSING";
        count=0
        for bundle in $bundles;
            do b=${bundle%"_to_template.trk"};
            echo $b;
            scil_volume_math.py lower_threshold_eq ${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/voxel_density_mask_voxel-norm_${b}_to_template.nii.gz 1.0 ${b}_mask.nii.gz -f;

            erosion=$(echo $erosions | awk -v n=$((count+1)) '{print $n}');
            echo "Erosion: $erosion";
            scil_volume_math.py erosion ${b}_mask.nii.gz ${erosion} ${b}_mask_eroded.nii.gz -f;
            ((count++));

        done;
    fi
done;
