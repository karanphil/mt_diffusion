#!/bin/bash

# The script assumes that preprocessing_dwi_pipeline.sh,
# preprocessing_t1_pipeline.sh, processing_fodf.sh, processing_tractogram.sh,
# processing_rbx_pipeline.sh, processing_bundles_pipeline.sh and
# processing_fixel_mtr_pipeline.sh have already been run.

# The first argument of the script is the target directory (full path)
target_dir=$1; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/full_processing"
# The second argument is the directory where the register_flow code lives (full path)
register_code_dir=$2; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/code/register_flow"
# The third argument is the register_flow directory to create and use (full path)
register_data_dir=$3; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/register_flow"
# The fourth argument is the T1 template file (has to be bet) (full path)
template_path=$4; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/mni_atlas/t1_template_bet.nii.gz"
# The fifth argument is the path to the singularity container (full path)
singularity_path=$5; # ex: "/home/local/USHERBROOKE/karp2601/data/stockage/mt-diff-mcgill/containers/scilus_2.2.0.sif"
cd $target_dir;
subs=$(ls -d hc*);

# ----------------------------Registration PROCESSING-------------------------
echo "Register-flow PROCESSING";
mkdir -p ${register_data_dir};
ln -s $template_path ${register_data_dir}/template_t1.nii.gz;
template_path="${register_data_dir}/template_t1.nii.gz";
for sub in $subs; 
    do echo $sub;

    mkdir -p ${register_data_dir}/input/${sub}/metrics;
    mkdir -p ${register_data_dir}/input/${sub}/tractograms;

    cd ${register_data_dir}/input/${sub};
    ln -sf ${target_dir}/${sub}/preprocessing_t1/register_natif/outputWarped.nii.gz t1.nii.gz

    cd ${register_data_dir}/input/${sub}/metrics;
    ln -s ${target_dir}/${sub}/fixel_analysis/voxel_density_mask_voxel-norm_*.nii.gz .;
    ln -s ${target_dir}/${sub}/fixel_mtr/fixel_mtr*.nii.gz .;
    ln -s ${target_dir}/${sub}/mtr/mtr*.nii.gz .;
    ln -s ${target_dir}/${sub}/afd_fixel/afd_fixel_*.nii.gz .;
    ln -s ${target_dir}/${sub}/fodf_metrics_mt_off/nufo.nii.gz .;
    ln -s ${target_dir}/${sub}/mtr/powder_averaged_mtr.nii.gz .;
    ln -s ${target_dir}/${sub}/dti/fa.nii.gz .;
    ln -s ${target_dir}/${sub}/dti/md.nii.gz .;

    cd ${register_data_dir}/input/${sub}/tractograms;
    ln -s ${target_dir}/${sub}/bundles/*.trk .;

done;

mkdir -p ${register_data_dir}/output;
cd ${register_data_dir}/output;

nextflow run ${register_code_dir}/main.nf --input ../input --template $template_path -with-singularity $singularity_path --trk_remove_invalid -resume;
