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
map_thr=1.0;
afd_thr=0.3; # This can be changed if needed
min_nvox=100; # This can be changed if needed
min_nb_subjects=5; # This can be changed if needed

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
            scil_bundle_label_map.py ${register_data_dir}/output/results_registration/${sub}/Trks_into_template_space/${bundle} ${target_dir}/average_all_scans/centroid_${b}.trk tmp --nb_pts $nb_sections -f;
            cp tmp/labels_map.nii.gz labels_${b}.nii.gz;
        
        done;
    fi
    rm -r tmp;

    # -------------------------Bundles masks PROCESSING------------------------
    cd ${target_dir}/${sub};
    mkdir -p masks_bundles;
    cd ${target_dir}/${sub}/masks_bundles;
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

    # -------------------------Track-profiles PROCESSING-----------------------
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
            nufo="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/nufo.nii_to_template.gz";
            labels="${target_dir}/${sub}/labels/labels_${b}.nii.gz";
            bundle_mask="${target_dir}/${sub}/masks_bundles/${b}_mask_eroded.nii.gz";

            python ${code_dir}/python_scripts/compute_track_profiles.py $mtr $fixel_mtr $labels ${b} $afd_fixel $nufo . --in_bundle_map $bundle_mask --map_threshold $map_thr --afd_threshold $afd_thr --min_nvox $min_nvox -f; 
        
        done;
    fi

done;

# ------------------------Track-profile plots PROCESSING-----------------------
cd $target_dir;
subs=$(ls -d hc[0-9][0-9]);
subs_scan=$(ls -d hc[0-9][0-9] 2>/dev/null | grep -v hc17 | grep -v hc28 | grep -v hc30)
subs_rescan=$(ls -d hc[0-9][0-9]r);

mkdir -p subject_wise_analysis;
cd ${target_dir}/subject_wise_analysis;

for bundle in $bundles;
    do b=${bundle%"_to_template.trk"};
    echo $b;
    # All scans
    in_mtr_profiles=$(for sub in $subs; do echo "${target_dir}/${sub}/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles=$(for sub in $subs; do echo "${target_dir}/${sub}/fixel_mtr_profile_${b}.txt"; done );
    in_nufo_profiles=$(for sub in $subs; do echo "${target_dir}/${sub}/nufo_profile_${b}.txt"; done );
    in_afd_profiles=$(for sub in $subs; do echo "${target_dir}/${sub}/afd_profile_${b}.txt"; done );
    # Only scans with rescans
    in_mtr_profiles_scan=$(for sub in $subs_scan; do echo "${target_dir}/${sub}/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles_scan=$(for sub in $subs_scan; do echo "${target_dir}/${sub}/fixel_mtr_profile_${b}.txt"; done );
    # in_nufo_profiles_scan=$(for sub in $subs_scan; do echo "${target_dir}/${sub}/nufo_profile_${b}.txt"; done );
    # in_afd_profiles_scan=$(for sub in $subs_scan; do echo "${target_dir}/${sub}/afd_profile_${b}.txt"; done );
    # Only rescans
    in_mtr_profiles_rescan=$(for sub in $subs_rescan; do echo "${target_dir}/${sub}/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles_rescan=$(for sub in $subs_rescan; do echo "${target_dir}/${sub}/fixel_mtr_profile_${b}.txt"; done );
    # in_nufo_profiles_rescan=$(for sub in $subs_rescan; do echo "${target_dir}/${sub}/nufo_profile_${b}.txt"; done );
    # in_afd_profiles_rescan=$(for sub in $subs_rescan; do echo "${target_dir}/${sub}/afd_profile_${b}.txt"; done );

    python ${code_dir}/python_scripts/plot_track_profiles_from_txt.py ${b} . --in_mtr_profiles_all $in_mtr_profiles --in_fixel_mtr_profiles_all $in_fixel_mtr_profiles --in_afd_fixel_profiles_all $in_afd_profiles --in_nufo_profiles_all $in_nufo_profiles --in_mtr_profiles_scan $in_mtr_profiles_scan --in_fixel_mtr_profiles_scan $in_fixel_mtr_profiles_scan --in_mtr_profiles_rescan $in_mtr_profiles_rescan --in_fixel_mtr_profiles_rescan $in_fixel_mtr_profiles_rescan --nb_sections $nb_sections --min_nb_subjects $min_nb_subjects -f;
    
done;