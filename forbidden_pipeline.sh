#!/bin/bash

# IMPORTANT: LAUNCH THIS SCRIPT FROM A FOLDER WHERE YOU WANT TO CREATE THE PROJECT FOLDER STRUCTURE.
# AKA inside mt-diff-mcgill folder.

# # Step 1
# main_dir="mt-diff-mcgill"; # Put any name you like
# mkdir $main_dir;
# cd $main_dir;
main_dir=$(pwd); # To get full paths for the future

# Step 2
code_name="code"; # Put any name you like
# mkdir $code_name;
# cd $code_name;
# git clone git@github.com:karanphil/mt_diffusion.git;
# cd mt-diffusion;
# git checkout mcgill_project; # It is important to use the mcgill_project branch.
code_dir="${main_dir}/${code_name}/mt_diffusion"; # Important to set this right

# Step 3
cd $main_dir;
working_dir="full_processing"; # Put any name you like
mkdir $working_dir;

# Step 4
original_data_dir="mt-diff-10peeps"; # IMPORTANT: Replace this with the actual path to the data
bash ${code_dir}/rename_files.sh $original_data_dir ${main_dir}/${working_dir};

# Step 5
container_dir="containers"; # Put any name you like
singularity_path="${main_dir}/${container_dir}/scilus_2.2.0.sif";
singularity build $singularity_path docker://scilus/scilpy:2.2.0;

# Step 6
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_dwi_pipeline.sh ${main_dir}/${working_dir} ${code_dir};

# Step 7
atlas_dir="${main_dir}/mni_atlas"; 
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_t1_pipeline.sh ${main_dir}/${working_dir} ${main_dir}/${atlas_dir};

# Step 8
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fodf_pipeline.sh ${main_dir}/${working_dir};

# Step 9
use_gpu="true"; # Set this to "false" if no gpu is available.
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_tractogram_pipeline.sh ${main_dir}/${working_dir} $use_gpu;

# Step 10
cd ${main_dir}/${code_name};
git clone git@github.com:scilus/rbx_flow.git;
rbx_code_dir="${main_dir}/${code_name}/rbx_flow";

# Step 11
rbx_data_dir="${main_dir}/rbx_flow"; # Put any name you like
rbx_atlas_dir="${main_dir}/rbx_atlas"; # Put the right path to the atlas
cd ${main_dir};
bash ${code_dir}/processing_rbx_pipeline.sh ${main_dir}/${working_dir} $rbx_code_dir $rbx_data_dir $rbx_atlas_dir $singularity_path;

# Step 12
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_bundles_pipeline.sh ${main_dir}/${working_dir} ${code_dir} $rbx_data_dir;

# Step 13
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fixel_mtr_pipeline.sh ${main_dir}/${working_dir} ${code_dir};

