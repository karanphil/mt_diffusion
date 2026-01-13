So the subject 25 (hc25 and hc25r) has very a thin CC, with hyperintensities on the T1 also. Thus, fsl-fast has trouble segmenting a proper WM mask in the middle part of the CC. To fix this, I computed a supplementary mask with the FA, using a threshold of FA>=0.5. I did the union of both wm_mask.nii.gz and this fa_thr_0p5.nii.gz. Here are the command:

!!! Don't forget to remove all files that you want to overwrite using the scripts, as we look for already existing data before launch !!!

```
# In the main directory of the project (mt-diff-mcgill), do:
main_dir=$(pwd -P);
working_dir="full_processing";
gpu_option="--use_gpu";
target_dir=${main_dir}/${working_dir};
container_dir="containers"; # Put any name you like
singularity_path="${main_dir}/${container_dir}/scilus_2.1.2.sif";
# Select one of hc25 or hc25r
sub=hc25;
# sub=hc25r;
# Go to the ${target_dir}/${sub}/preprocessing_t1/register_natif folder.
cd ${target_dir}/${sub}/preprocessing_t1/register_natif;
scil_volume_math.py lower_threshold_eq ${target_dir}/${sub}/dti/fa.nii.gz 0.5 fa_thr_0p5.nii.gz;
scil_volume_math.py union wm_mask.nii.gz fa_thr_0p5.nii.gz wm_mask_for_tractography.nii.gz
```

Then, I ran the tractography using this mask.

```
# Go to the ${target_dir}/${sub}/tractography folder.
cd ${target_dir}/${sub}/tractography;
scil_volume_resample.py ${target_dir}/${sub}/preprocessing_t1/register_natif/wm_mask_for_tractography.nii.gz ${target_dir}/${sub}/dwi_for_tractography/wm_mask_extended.nii.gz --voxel_size 1 --interp nn -f;
fodf="${target_dir}/${sub}/fodf_for_tractography/fodf.nii.gz";
wm_mask="${target_dir}/${sub}/dwi_for_tractography/wm_mask_extended.nii.gz";
scil_tracking_local.py $fodf $wm_mask $wm_mask local_tracking.trk $gpu_option --npv 10 -f;
scil_tractogram_remove_invalid.py local_tracking.trk local_tracking.trk --remove_single_point -f;
```

Now, run SIFT2:

```
fodf_tournier="${target_dir}/${sub}/fodf_for_tractography/fodf_tournier.nii.gz";
scil_sh_convert.py $fodf $fodf_tournier descoteaux07_legacy tournier07 -f;
tractogram_tck="${target_dir}/${sub}/tractography/local_tracking.tck";
scil_tractogram_convert.py $tractogram $tractogram_tck -f;
tcksift2 $tractogram_tck $fodf_tournier sift2_weights.txt -force;
scil_tractogram_dps_math.py $tractogram import "sift2" --in_dps_file sift2_weights.txt --out_tractogram $tractogram -f;
rm $tractogram_tck $fodf_tournier;
```

If rbx_flow was already ran, update the tractograms in rbx_flow, and resume it:

```
rbx_flow_version="1.3.0";
rbx_data_dir="${main_dir}/rbx_flow";
rbx_atlas_dir="${main_dir}/rbx_atlas";
sub=hc25;
tractogram="${target_dir}/${sub}/tractography/local_tracking.trk";
ln -sf $tractogram ${rbx_data_dir}/input/${sub}/local_tracking.trk;
sub=hc25r;
tractogram="${target_dir}/${sub}/tractography/local_tracking.trk";
ln -sf $tractogram ${rbx_data_dir}/input/${sub}/local_tracking.trk;

cd ${rbx_data_dir}/output;
nextflow run scilus/rbx_flow -r $rbx_flow_versions --input ../input --atlas_directory $rbx_atlas_dir -with-singularity $singularity_path --register_processes 8 --rbx_processes 8 -resume;
```