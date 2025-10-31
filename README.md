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
code_dir="code"; # Put any name you like
mkdir $code_dir;
cd $code_dir;
git clone git@github.com:karanphil/mt_diffusion.git;
cd mt-diffusion;
git checkout mcgill_project; # It is important to use the mcgill_project branch.
code_dir="${code_dir}/mt_diffusion"; # Important to set this right
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
bash ${main_dir}/${code_dir}/rename_files.sh $original_data_dir ${main_dir}/${working_dir};
```

Now, compute the pre-processing of DWI data:

```
bash ${main_dir}/${code_dir}/preprocessing_dwi_pipeline.sh ${main_dir}/${working_dir} ${main_dir}/${code_dir};
```

When this is finished (after many hours), compute the pre-processing of T1 data. You will need a T1 atlas (named t1_template.nii.gz) with a brain probability map (named t1_brain_probability_map.nii.gz). Put them in the $main_dir directory in a folder named "mni_atlas". Then:

```
atlas_dir="${main_dir}/mni_atlas"; 
bash ${main_dir}/${code_dir}/preprocessing_t1_pipeline.sh ${main_dir}/${working_dir} ${main_dir}/${atlas_dir};
```