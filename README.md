These are the instructions for running the whole MT-diffusion pipeline.

First, you should create a repository to host the project: 

```
main_dir="mt-diff-mcgill"; # Put any name you like
mkdir $main_dir;
cd $main_dir;
main_dir=$(pwd); # To get full paths for the future
```

Next, create a folder for your code and clone the mt-diffusion code repository in it:

```
code_name="code"; # Put any name you like
mkdir $code_name;
cd $code_name;
git clone git@github.com:karanphil/mt_diffusion.git;
cd mt-diffusion;
git checkout mcgill_project; # It is important to use the mcgill_project branch.
code_dir="${main_dir}/${code_name}/mt_diffusion"; # Important to set this right
```

Then, go back to the $main_dir directory and create the folder where all the data will be stored:

```
cd $main_dir;
working_dir="full_processing"; # Put any name you like
mkdir $working_dir;
```

Then, it is time to copy the data the $working_dir directory, using the rename_files.sh script:

```
original_data_dir="mt-diff-10peeps"; # IMPORTANT: Replace this with the actual path to the data
bash ${code_dir}/rename_files.sh $original_data_dir ${main_dir}/${working_dir};
```

For the next steps, it is suggested to use a singularity. It can be downloaded from Docker Hub (https://hub.docker.com/r/scilus/scilus), as the scilus 2.2.0, or simply using the follow command lines:

```
container_dir="containers"; # Put any name you like
singularity_path="${main_dir}/${container_dir}/scilus_2.2.0.sif";
singularity build $singularity_path docker://scilus/scilpy:2.2.0;
```

Now, compute the pre-processing of DWI data:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_dwi_pipeline.sh ${main_dir}/${working_dir} ${code_dir};
```

When this is finished (after many hours), compute the pre-processing of T1 data. You will need a T1 atlas (named t1_template.nii.gz) with a brain probability map (named t1_brain_probability_map.nii.gz). Put them in the $main_dir directory in a folder named "mni_atlas". Then:

```
atlas_dir="${main_dir}/mni_atlas"; 
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_t1_pipeline.sh ${main_dir}/${working_dir} ${main_dir}/${atlas_dir};
```

After that, compute the fODFs with:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fodf_pipeline.sh ${main_dir}/${working_dir};
```

Then, compute the tractogram with:

```
use_gpu="true"; # Set this to "false" if no gpu is available.
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_tractogram_pipeline.sh ${main_dir}/${working_dir} $use_gpu;
```

Then, run rbx_flow. You will need to first clone the rbx_flow repository:

```
cd ${main_dir}/${code_name};
git clone git@github.com:scilus/rbx_flow.git;
rbx_code_dir="${main_dir}/${code_name}/rbx_flow";
```

Second, a bundles atlas is needed. Add the path to this atlas and run the rbx pipeline:

```
rbx_data_dir="${main_dir}/rbx_flow"; # Put any name you like
rbx_atlas_dir="${main_dir}/rbx_atlas"; # Put the right path to the atlas
cd ${main_dir};
bash ${code_dir}/processing_rbx_pipeline.sh ${main_dir}/${working_dir} $rbx_code_dir $rbx_data_dir $rbx_atlas_dir $singularity_path;
```

After computing the bundles with rbx_flow, it is time to bring the SIFT2 weights back from the tractogram to the bundles, and then run a bundle fixel analysis (fixel density and other stuff):

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_bundles_pipeline.sh ${main_dir}/${working_dir} ${code_dir} $rbx_data_dir;
```