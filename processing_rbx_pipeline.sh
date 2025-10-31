#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh and
# processing_tractogram.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the directory where the rbx_flow code lives (full path)
rbx_code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/rbx_flow"
# The third argument is the rbx_flow directory to create and use (full path)
rbx_data_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/rbx_flow"
# The fourth argument is the bundles atlas directory (full path)
rbx_atlas_dir=$4; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/rbx_atlas"
# The fifth argument is the path to the singularity container (full path)
singularity_path=$5; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/containers/scilus_2.2.0.sif"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------Bundles PROCESSING-------------------------
echo "RBX-Flow PROCESSING";
for sub in $subs; 
    do echo $sub;

    fa="${target_dir}/${sub}/dwi_for_tractography/fa.nii.gz";
    tractogram="${target_dir}/${sub}/tractography/local_tracking.trk";

    mkdir -p ${rbx_data_dir}/input/${sub};
    ln -s $fa ${rbx_data_dir}/input/${sub}/fa.nii.gz;
    ln -s $tractogram ${rbx_data_dir}/input/${sub}/local_tracking.trk;

done;

mkdir -p ${rbx_data_dir}/output;
cd ${rbx_data_dir}/output;

nextflow run ${rbx_code_dir}/main.nf --input ../input --atlas_directory $rbx_atlas_dir -with-singularity $singularity_path --register_processes 8 --rbx_processes 8 -resume;
