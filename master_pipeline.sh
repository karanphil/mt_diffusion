#!/bin/bash

# IMPORTANT: LAUNCH THIS SCRIPT FROM A FOLDER WHERE YOU WANT TO CREATE THE PROJECT FOLDER STRUCTURE.
# AKA inside mt-diff-mcgill folder.

# VERSIONS
rbx_flow_version="1.3.0";
register_flow_version="0.2.0";

# # Step 1
# main_dir="mt-diff-mcgill"; # Put any name you like
# mkdir $main_dir;
# cd $main_dir;
main_dir=$(pwd -P); # To get full paths for the future

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
mkdir -p $working_dir;

# Step 4
original_data_dir="mt-diff-10peeps"; # IMPORTANT: Replace this with the actual path to the data
bash ${code_dir}/rename_files.sh $original_data_dir ${main_dir}/${working_dir};

# Step 5
container_dir="containers"; # Put any name you like
singularity_path="${main_dir}/${container_dir}/scilus_2.1.2.sif";
singularity build $singularity_path docker://scilus/scilus:2.1.2;

# Step 6
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_dwi_pipeline.sh ${main_dir}/${working_dir} ${code_dir};

# Step 7
atlas_dir="${main_dir}/mni_atlas"; 
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_t1_pipeline.sh ${main_dir}/${working_dir} ${main_dir}/${atlas_dir};

# Step 8
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fodf_pipeline.sh ${main_dir}/${working_dir};

# Step 9
use_gpu="true"; # Set this to "false" if no gpu is available.
nvidia_opt="--nv"; # Comment if no gpu is available.
singularity exec $nvidia_opt -B $main_dir $singularity_path bash ${code_dir}/processing_tractogram_pipeline.sh ${main_dir}/${working_dir} $use_gpu;

# Step 10
rbx_data_dir="${main_dir}/rbx_flow"; # Put any name you like
rbx_atlas_dir="${main_dir}/rbx_atlas"; # Put the right path to the atlas
cd ${main_dir};
bash ${code_dir}/processing_rbx_pipeline.sh ${main_dir}/${working_dir} $rbx_flow_version $rbx_data_dir $rbx_atlas_dir $singularity_path;

# Step 11
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_bundles_pipeline.sh ${main_dir}/${working_dir} ${code_dir} $rbx_data_dir;

# Step 12
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fixel_mtr_pipeline.sh ${main_dir}/${working_dir} ${code_dir};

# Step 13
register_data_dir="${main_dir}/register_flow";
template_path="${main_dir}/mni_atlas/t1_template_bet.nii.gz"; # Put the right path to the template
cd ${main_dir};
bash ${code_dir}/processing_registration_pipeline.sh ${main_dir}/${working_dir} $register_flow_version $register_data_dir $template_path $singularity_path;

# Step 14
register_data_dir="${main_dir}/register_flow";
cd ${main_dir};
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_averages_pipeline.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};

# Step 15
average_dir=${main_dir}/${working_dir}/average_all_scans;
cd ${main_dir};
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_average_pipeline.sh ${average_dir} ${code_dir};

# Step 16
register_data_dir="${main_dir}/register_flow";
cd ${main_dir};
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_subwise_pipeline.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};

# Step 17
register_data_dir="${main_dir}/register_flow";
cd ${main_dir};
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_bundles_crossing_stats.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};

# Step 18
register_data_dir="${main_dir}/register_flow";
cd ${main_dir};
# singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_stats.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
# The singularity's version of scipy does not have the nan_policy option for the shapiro test. Please run this without singularity. If the scilpy dependencies are problematic, remove those functions from the script or install scilpy.
bash ${code_dir}/processing_track_profiles_stats.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
