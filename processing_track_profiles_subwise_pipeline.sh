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
# The fourth argument is wether to include overlap and significance info in the plots
include_extra_info=$4; # ex: True or False
cd $target_dir;
subs=$(ls -d hc*);

nb_sections=20; # This can be changed if needed
map_thr=1.0; # Take all the eroded masks
afd_thr=0.3; # This can be changed if needed
min_nvox=100; # This can be changed if needed
min_nb_subjects=5; # This can be changed if needed

overlap_cmd="";
significance_cmd="";
if [ $include_extra_info == "True" ]; then
    overlap_cmd="--in_overlap_txt ${register_data_dir}/output/processing_registration/crossing_bundles_labels_important.txt";
    significance_cmd="--in_significance_txt ${register_data_dir}/output/processing_registration/significant_track_profiles_crossings.txt";
fi

cd ${register_data_dir}/output/results_registration/hc18/Trks_into_template_space;
bundles=$(ls *.trk);
for sub in $subs; 
    do echo $sub;
    # -------------------------Track-profiles PROCESSING-----------------------
    cd ${register_data_dir}/output;
    mkdir -p processing_registration/${sub}/track_profiles;
    cd ${register_data_dir}/output/processing_registration/${sub}/track_profiles;
    if [ ! -f "fixel_mtr_profile_AF_L.txt" ]; then
        echo "Track-profiles PROCESSING";
        for bundle in $bundles;
            do b=${bundle%"_to_template.trk"};
            echo $b;

            mtr="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/mtr_${b}_to_template.nii.gz";
            fixel_mtr="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/fixel_mtr_${b}_to_template.nii.gz";
            afd_fixel="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/afd_fixel_${b}_to_template.nii.gz";
            nufo="${register_data_dir}/output/results_registration/${sub}/Metrics_into_template_space/nufo_to_template.nii.gz";
            labels="${register_data_dir}/output/processing_registration/${sub}/labels/labels_${b}.nii.gz";
            # bundle_mask="${register_data_dir}/output/processing_registration/${sub}/bundles_masks/${b}_mask_eroded.nii.gz";
            bundle_mask="${register_data_dir}/output/processing_registration/${sub}/bundles_masks/${b}_mask.nii.gz";

            python ${code_dir}/python_scripts/compute_track_profiles.py $mtr $fixel_mtr $labels ${b} $afd_fixel $nufo . --in_bundle_map $bundle_mask --map_threshold $map_thr --afd_threshold $afd_thr --min_nvox $min_nvox --nb_sections $nb_sections -f; 
        
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

echo "Track-profiles PLOTTING";
for bundle in $bundles;
    do b=${bundle%"_to_template.trk"};
    echo $b;
    # All scans
    in_mtr_profiles=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/fixel_mtr_profile_${b}.txt"; done );
    in_nufo_profiles=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/nufo_profile_${b}.txt"; done );
    in_afd_profiles=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/afd_profile_${b}.txt"; done );
    # Only scans with rescans
    in_mtr_profiles_scan=$(for sub in $subs_scan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles_scan=$(for sub in $subs_scan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/fixel_mtr_profile_${b}.txt"; done );
    # in_nufo_profiles_scan=$(for sub in $subs_scan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/nufo_profile_${b}.txt"; done );
    # in_afd_profiles_scan=$(for sub in $subs_scan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/afd_profile_${b}.txt"; done );
    # Only rescans
    in_mtr_profiles_rescan=$(for sub in $subs_rescan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles_rescan=$(for sub in $subs_rescan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/fixel_mtr_profile_${b}.txt"; done );
    # in_nufo_profiles_rescan=$(for sub in $subs_rescan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/nufo_profile_${b}.txt"; done );
    # in_afd_profiles_rescan=$(for sub in $subs_rescan; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/afd_profile_${b}.txt"; done );

    python ${code_dir}/python_scripts/plot_track_profiles_from_subjects.py ${b} . --in_mtr_profiles_all $in_mtr_profiles --in_fixel_mtr_profiles_all $in_fixel_mtr_profiles --in_afd_fixel_profiles_all $in_afd_profiles --in_nufo_profiles_all $in_nufo_profiles --in_mtr_profiles_scan $in_mtr_profiles_scan --in_fixel_mtr_profiles_scan $in_fixel_mtr_profiles_scan --in_mtr_profiles_rescan $in_mtr_profiles_rescan --in_fixel_mtr_profiles_rescan $in_fixel_mtr_profiles_rescan --nb_sections $nb_sections --min_nb_subjects $min_nb_subjects $overlap_cmd $significance_cmd -f;
    
done;