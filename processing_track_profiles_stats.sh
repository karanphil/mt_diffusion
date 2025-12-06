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

cd ${register_data_dir}/output/results_registration/hc18/Trks_into_template_space;
bundles=$(ls *.trk);
cd $target_dir;
subs=$(ls -d hc[0-9][0-9]);
cd ${register_data_dir}/output/processing_registration/;
rm ${register_data_dir}/output/processing_registration/significant_track_profiles_crossings.txt; # Remove previous file if exists

for bundle in $bundles;
    do b=${bundle%"_to_template.trk"};
    echo $b;
    # All scans
    in_mtr_profiles=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/mtr_profile_${b}.txt"; done );
    in_fixel_mtr_profiles=$(for sub in $subs; do echo "${register_data_dir}/output/processing_registration/${sub}/track_profiles/fixel_mtr_profile_${b}.txt"; done );

    python ${code_dir}/python_scripts/compute_track_profiles_stats.py ${b} ${register_data_dir}/output/processing_registration/significant_track_profiles_crossings.txt --in_mtr_profiles $in_mtr_profiles --in_fixel_mtr_profiles $in_fixel_mtr_profiles -f;
    
done;