#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh,
# processing_rbx_pipeline.sh, processing_bundles_pipeline.sh,
# processing_fixel_mtr_pipeline.sh and processing_registration_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"
# The third argument is the register_flow directory (full path)
register_data_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/register_flow"
cd $target_dir;
# subs=$(ls -d hc*);
subs=$(ls -d hc[0-9][0-9]);

# If we ever need to compare the population average scan-rescan, do the code twice with the following subset of subjects.
# subs_scan=$(ls -d hc[0-9][0-9] 2>/dev/null | grep -v hc17 | grep -v hc28 | grep -v hc30) # Command for removing subjects without rescan
# subs_rescan=$(ls -d hc[0-9][0-9]r);

nb_sections=20; # This can be changed if needed

# ----------------------------Average images PROCESSING------------------------
cd ${register_data_dir}/output/results_registration/hc18/Metrics_into_template_space;
images=$(ls);
echo "Average images PROCESSING";
cd $target_dir;
mkdir -p average;
cd ${target_dir}/average;
for image in $images; 
    do echo $image;

    inputs=$(for sub in $subs; do echo "${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/${image}"; done );
    image_name=${image%"_to_template.nii.gz"}

    if [[ "$image" == *"voxel_density_mask"* ]]; then
        mean_median_option="";
    else
        mean_median_option="--median";
    fi

    python ${code_dir}/python_scripts/compute_average_across_acquisitions.py $inputs ${image_name}.nii.gz $mean_median_option -f;

done;

# ----------------------------Average bundles PROCESSING-----------------------
cd ${register_data_dir}/output/results_registration/hc18/Trks_into_template_space;
bundles=$(ls *.trk);
echo "Average bundles PROCESSING";
mkdir -p tmp;
for bundle in $bundles;
    do b=${bundle%"_to_template.trk"};
    echo $b;

    inputs=$(for sub in $subs; do echo "${register_data_dir}/output/results_registration/${sub}/Trks_into_template_space/${bundle}"; done );
    scil_tractogram_math union $inputs united_${b}.trk -f;

    scil_tractogram_compress united_${b}.trk united_${b}_compressed.trk -f;
    scil_bundle_alter_to_target_dice united_${b}_compressed.trk united_${b}_altered.trk --subsample --min_dice 0.95 -f;
    mv united_${b}_altered.trk united_${b}.trk;
    rm united_${b}_compressed.trk;
    scil_bundle_uniformize_endpoints united_${b}.trk united_${b}_tmp.trk --auto -f;
    mv united_${b}_tmp.trk united_${b}.trk;
    scil_bundle_compute_centroid united_${b}.trk centroid_${b}.trk --nb_points $nb_sections -f;
    scil_bundle_label_map united_${b}.trk centroid_${b}.trk tmp --nb_pts $nb_sections -f;
    cp tmp/labels_map.nii.gz labels_${b}.nii.gz;

done;
rm -r tmp;
