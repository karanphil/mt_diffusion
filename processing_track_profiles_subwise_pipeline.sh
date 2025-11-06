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
subs=$(ls -d hc*);
nb_sections=20; # This can be changed if needed

cd ${register_data_dir}/output/results_registration/hc18/Trks_into_template_space;
bundles=$(ls *.trk);
for sub in $subs; 
    do echo $sub;
    # ----------------------------Labels PROCESSING-------------------------
    cd ${target_dir}/${sub};
    mkdir -p labels;
    cd ${target_dir}/${sub}/labels;
    mkdir -p tmp;
    if [ ! -f "labels_AF_L.nii.gz" ]; then
        echo "Labels PROCESSING";
        for bundle in $bundles;
            do b=${bundle%"_to_template.trk"};
            echo $b;
            scil_bundle_label_map ${register_data_dir}/output/results_registration/${sub}/Trks_into_template_space/${bundle} ${target_dir}/average/centroid_${b}.trk tmp --nb_pts $nb_sections -f;
            cp tmp/labels_map.nii.gz labels_${b}.nii.gz;
        
        done;
    fi
    rm -r tmp;

    # ----------------------------Track-profiles PROCESSING-------------------------
    cd ${target_dir}/${sub};
    mkdir -p track_profiles;
    cd ${target_dir}/${sub}/track_profiles;
    if [ ! -f "fixel_mtr_profile_AF_L.nii.gz" ]; then
        echo "Track-profiles PROCESSING";
        for bundle in $bundles;
            do b=${bundle%"_to_template.trk"};
            echo $b;

            mtr="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/mtr_${b}_to_template.nii.gz";
            fixel_mtr="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/fixel_mtr_${b}_to_template.nii.gz";
            afd_fixel="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/afd_fixel_${b}_to_template.nii.gz";
            nufo="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/nufo_to_template.nii.gz";
            bundle_map="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/voxel_density_mask_voxel-norm_${b}_to_template.nii.gz";

            python ${code_dir}/python_scripts/compute_track_profiles.py $mtr $fixel_mtr ../labels/labels_${b}.nii.gz $afd_fixel $nufo . --in_bundle_map $bundle_map --bundle_name ${b} -f; 
        
        done;
    fi

done;

# ------------------------Track-profile plots PROCESSING-----------------------
cd $target_dir;
subs_scan=$(ls -d hc[0-9][0-9]);

mkdir -p

for bundle in $bundles;
    do b=${bundle%"_to_template.trk"};
    echo $b;
    in_mtr_profiles=$(for d in $subs_scan; do echo "$d/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles=$(for d in $subs_scan; do echo "$d/fixel_mtr_profile_${b}.txt"; done );
    in_nufo_profiles=$(for d in $subs_scan; do echo "$d/nufo_profile_${b}.txt"; done );
    in_afd_profiles=$(for d in $subs_scan; do echo "$d/afd_profile_${b}.txt"; done );
    python ${code_dir}/python_scripts/plot_track_profiles_from_txt.py . --in_mtr_profiles $in_mtr_profiles --in_fixel_mtr_profiles $in_fixel_mtr_profiles --in_nufo_profiles $in_nufo_profiles --in_afd_profiles $in_afd_profiles --bundle_name ${b} --variance -f;
    
done;