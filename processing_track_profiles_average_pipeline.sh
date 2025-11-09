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
# The third argument is the average directory (full path)
average_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/average"

nb_sections=20; # This can be changed if needed

cd ${average_dir};
bundles=$(ls *.trk);

# ------------------------Track-profile plots PROCESSING-----------------------
cd $target_dir;
mkdir -p average_analysis;
cd ${target_dir}/average_analysis;

for bundle in $bundles;
    do b=${bundle#"united_"};
    b=${b%".trk"};
    echo $b;
    python ${code_dir}/python_scripts/plot_track_profiles.py ${average_dir}/mtr_${b}.nii.gz ${average_dir}/fixel_mtr_${b}.nii.gz ${average_dir}/labels_${b}.nii.gz ${average_dir}/afd_fixel_${b}.nii.gz ${average_dir}/nufo.nii.gz . --in_bundle_map ${average_dir}/voxel_density_mask_voxel-norm_${b}_to_template.nii.gz --bundle_name ${b} --map_threshold 1.0 --afd_threshold 0.3;
    
done;