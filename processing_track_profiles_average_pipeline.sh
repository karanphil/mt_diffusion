#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh,
# processing_rbx_pipeline.sh, processing_bundles_pipeline.sh,
# processing_fixel_mtr_pipeline.sh, processing_registration_pipeline.sh and
# processing_averages_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/average_all_scans"
# The second argument is the code directory (full path).
code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/mt-diffusion"

nb_sections=20; # This can be changed if needed
map_thr=1.0; # This can be changed if needed
afd_thr=0.3; # This can be changed if needed
min_nvox=100; # This can be changed if needed

# ------------------------Track-profile plots PROCESSING-----------------------
cd $target_dir;
bundles=$(ls united_*.trk);
mkdir -p track_profiles;
cd ${target_dir}/track_profiles;

for bundle in $bundles;
    do b=${bundle#"united_"};
    b=${b%".trk"};
    echo $b;
    python ${code_dir}/python_scripts/compute_track_profiles.py ${target_dir}/mtr_${b}.nii.gz ${target_dir}/fixel_mtr_${b}.nii.gz ${target_dir}/labels_${b}.nii.gz ${b} ${target_dir}/afd_fixel_${b}.nii.gz ${target_dir}/nufo.nii.gz . --nb_sections $nb_sections --map_threshold $map_thr --afd_threshold $afd_thr --min_nvox $min_nvox;

    python ${code_dir}/python_scripts/plot_track_profiles_from_average.py mtr_profile_${b}.txt fixel_mtr_profile_${b}.txt ${b} nufo_profile.txt afd_profile_${b}.txt . --nb_sections $nb_sections;
    
done;