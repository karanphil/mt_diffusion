These are the instructions for running the whole MT-diffusion pipeline.

Important note: If you ever want to comment out a section of a pipeline, always comment stuff inside IF statements, never outside.

First, set the versions for rbx_flow and register_flow:

```
rbx_flow_version="1.3.0";
register_flow_version="0.2.0";
```

Then, you should create a repository to host the project: 

```
main_dir="mt-diff-mcgill"; # Put any name you like
mkdir $main_dir;
cd $main_dir;
main_dir=$(pwd -P); # To get full paths for the future
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
singularity_path="${main_dir}/${container_dir}/scilus_2.1.2.sif";
singularity build $singularity_path docker://scilus/scilus:2.1.2;
```

Now, compute the pre-processing of DWI data:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_dwi_pipeline.sh ${main_dir}/${working_dir} ${code_dir};
```

When this is finished (after many hours), compute the pre-processing of T1 data. You will need a T1 atlas (named t1_template.nii.gz) with a brain probability map (named t1_brain_probability_map.nii.gz). Put them in the $main_dir directory in a folder named "mni_atlas". Then:

```
atlas_dir="${main_dir}/mni_atlas"; 
singularity exec -B $main_dir $singularity_path bash ${code_dir}/preprocessing_t1_pipeline.sh ${main_dir}/${working_dir} ${atlas_dir};
```

After that, compute the fODFs with:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fodf_pipeline.sh ${main_dir}/${working_dir};
```

Then, compute the tractogram with:

```
use_gpu="true"; # Set this to "false" if no gpu is available.
nvidia_opt="--nv"; # Comment if no gpu is available.
singularity exec $nvidia_opt -B $main_dir $singularity_path bash ${code_dir}/processing_tractogram_pipeline.sh ${main_dir}/${working_dir} $use_gpu;
```

Then, run rbx_flow. A bundles atlas is needed. Add the path to this atlas and run the rbx pipeline:

```
rbx_data_dir="${main_dir}/rbx_flow";
rbx_atlas_dir="${main_dir}/rbx_atlas"; # Put the right path to the atlas
cd ${main_dir};
bash ${code_dir}/processing_rbx_pipeline.sh ${main_dir}/${working_dir} $rbx_flow_version $rbx_data_dir $rbx_atlas_dir $singularity_path;
```

After computing the bundles with rbx_flow, it is time to bring the SIFT2 weights back from the tractogram to the bundles, and then run a bundle fixel analysis (fixel density and other stuff):

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_bundles_pipeline.sh ${main_dir}/${working_dir} ${code_dir} $rbx_data_dir;
```

Next, compute the fixel-MTR, PA-MTR and project them to the bundles:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_fixel_mtr_pipeline.sh ${main_dir}/${working_dir} ${code_dir};
```

Now that everything is computed per subject, it is time to register them all to the MNI template:

```
register_data_dir="${main_dir}/register_flow";
template_path="${main_dir}/mni_atlas/t1_template_bet.nii.gz"; # Put the right path to the template
cd ${main_dir};
bash ${code_dir}/processing_registration_pipeline.sh ${main_dir}/${working_dir} $register_flow_version $register_data_dir $template_path $singularity_path;
```

Next, compute the population average of all images and bundles:

```
register_data_dir="${main_dir}/register_flow";
cd ${main_dir};
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_averages_pipeline.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```

Now, if you want, compute and plot the track-profiles for the average data (this was removed from the final analysis):

```
average_dir=${main_dir}/${working_dir}/average_all_scans;
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_average_pipeline.sh ${average_dir} ${code_dir};
```

Next, compute the labels and bundle masks for each subjects using the average centroids:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_labels_per_subjects_pipeline.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```

Then, compute and plot the track-profiles for the between/within-subject analysis:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_subwise_pipeline.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```

Now, a matrix of percentage overlap comparing the overlap of each section of each bundle can be computed if wanted:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_bundles_crossing_stats.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```

Furthermore, the significance of the differences between MTR and fixel-MTR in track-profiles can be evaluated:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_stats.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```

If the Shapiro test for normality is uncommented from the python script, it will not run with the singularity. Simply install the necessary packages and run instead:

```
bash ${code_dir}/processing_track_profiles_stats.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```

If you want to add the important crossing regions and significant sections to the track-profiles plots, rerun this command:

```
singularity exec -B $main_dir $singularity_path bash ${code_dir}/processing_track_profiles_subwise_pipeline.sh ${main_dir}/${working_dir} ${code_dir} ${register_data_dir};
```
