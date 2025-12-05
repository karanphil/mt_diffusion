#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh,
# processing_rbx_pipeline.sh, processing_bundles_pipeline.sh,
# processing_fixel_mtr_pipeline.sh, processing_registration_pipeline.sh,
# processing_averages_pipeline.sh,
# processing_track_profiles_average_pipeline.sh and
# processing_track_profiles_subwise_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"
# The third argument is the register_flow directory (full path)
register_data_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/register_flow"
cd $target_dir;
subs=$(ls -d hc*);

nb_sections=20; # This can be changed if needed
map_thr=1.0; # Take all the eroded masks
afd_thr=0.3; # This can be changed if needed
min_nvox=100; # This can be changed if needed
highlight_threshold=0.2; # This can be changed if needed

in_bundles="AF_L AF_R CC_1 CC_2a CC_2b CC_3 CC_4 CC_5 CC_6 CC_7 CG_L CG_R CST_L CST_R ILF_L ILF_R MCP OR_L OR_R SLF_1_L SLF_1_R SLF_2_L SLF_2_R SLF_3_L SLF_3_R"; # Removing some bundles to simplify the plots

# ------------------------Crossing bundles labels PROCESSING-------------------
echo "Crossing bundles labels PROCESSING";
cd ${register_data_dir}/output;
for sub in $subs; 
    do echo $sub;
        if [ ! -f "processing_registration/${sub}/crossing_bundles_labels.json" ]; then
            python ${code_dir}/python_scripts/compute_crossing_bundles_labels.py processing_registration/${sub}/crossing_bundles_labels.json --in_labels processing_registration/${sub}/labels/labels_* --in_afd_fixel results_registration/${sub}/Metrics_into_template_space/afd_fixel_* --in_bundle_map results_registration/${sub}/Metrics_into_template_space/voxel_density_mask_voxel-norm_* --map_threshold $map_thr --afd_threshold $afd_thr --min_nvox $min_nvox --nb_sections $nb_sections -f;
        fi
done;

# ------------------------Track-profile plots PROCESSING-----------------------
cd $target_dir;
subs=$(ls -d hc[0-9][0-9]);
cd ${register_data_dir}/output/processing_registration;

echo "Crossing bundles labels PLOTTING";
in_jsons=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/crossing_bundles_labels.json"; done );
python ${code_dir}/python_scripts/plot_crossing_bundles_labels.py ${register_data_dir}/output/processing_registration/crossing_bundles_labels.png --in_jsons $in_jsons --bundles_names $in_bundles --nb_sections $nb_sections --highlight_threshold $highlight_threshold --save_txt ${register_data_dir}/output/processing_registration/crossing_bundles_labels.txt -f;
