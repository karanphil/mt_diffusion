mkdir acquisition1; cd acquisition1;
# denoising DWI with mrtrix LPCA (could be optional)
# you could also average your 3 acquisitions...
mkdir denoised_data; cd denoised_data;
dwidenoise ../../raw_data/Anthony_05_19_2022_first_test_WIP_MT-diff_20220519105946_401.nii dwi.nii.gz;
cp ../../raw_data/Anthony_05_19_2022_first_test_WIP_MT-diff_20220519105946_401.bval dwi.bval;
cp ../../raw_data/Anthony_05_19_2022_first_test_WIP_MT-diff_20220519105946_401.bvec dwi.bvec;
cd ../;

# Removing bad b0s
mkdir split_data; cd split_data;
# scil_split_image.py ../denoised_data/dwi.nii.gz ../denoised_data/dwi.bval ../denoised_data/dwi.bvec 3 4 36 39 40 --out_dwi dwi_b0_trash_mt_off.nii.gz dwi_b0_mt_off.nii.gz dwi_b1500_mt_off.nii.gz dwi_b0_trash_mt_on.nii.gz dwi_b0_mt_on.nii.gz dwi_b1500_mt_on.nii.gz --out_bval dwi_b0_trash_mt_off.bval dwi_b0_mt_off.bval dwi_b1500_mt_off.bval dwi_b0_trash_mt_on.bval dwi_b0_mt_on.bval dwi_b1500_mt_on.bval --out_bvec dwi_b0_trash_mt_off.bvec dwi_b0_mt_off.bvec dwi_b1500_mt_off.bvec dwi_b0_trash_mt_on.bvec dwi_b0_mt_on.bvec dwi_b1500_mt_on.bvec
scil_split_image.py ../denoised_data/dwi.nii.gz ../denoised_data/dwi.bval ../denoised_data/dwi.bvec dwi 3 4 36 39 40;
mkdir mt_off; cd mt_off;
scil_concatenate_dwi.py dwi.nii.gz dwi.bval dwi.bvec --in_dwis ../dwi_3_3.nii.gz ../dwi_4_35.nii.gz --in_bvals ../dwi_3_3.bval ../dwi_4_35.bval --in_bvecs ../dwi_3_3.bvec ../dwi_4_35.bvec -f;
cd ../;
mkdir mt_on; cd mt_on;
scil_concatenate_dwi.py dwi.nii.gz dwi.bval dwi.bvec --in_dwis ../dwi_39_39.nii.gz ../dwi_40_71.nii.gz --in_bvals ../dwi_39_39.bval ../dwi_40_71.bval --in_bvecs ../dwi_39_39.bvec ../dwi_40_71.bvec -f;
cd ../../;

# Eddy-topup
mkdir eddy_topup; cd eddy_topup
# MT OFF
mkdir mt_off; cd mt_off
# Extract b0 & brain extraction
scil_extract_b0.py ../../split_data/mt_off/dwi.nii.gz ../../split_data/mt_off/dwi.bval ../../split_data/mt_off/dwi.bvec dwi_mean_b0.nii.gz --b0_thr 20 --mean
bet dwi_mean_b0.nii.gz dwi_mean_b0_bet.nii.gz -m -R
# Register rev_b0 to mean_b0
antsRegistrationSyNQuick.sh -d 3 -f dwi_mean_b0.nii.gz -m ../../../raw_data/Anthony_05_19_2022_first_test_WIP_rev_b0_20220519105946_501.nii -o output -t r
mv outputWarped.nii.gz dwi_rev_b0_warped.nii.gz
# Topup
scil_prepare_topup_command.py dwi_mean_b0.nii.gz dwi_rev_b0_warped.nii.gz --out_script --out_prefix topup_results --out_b0s dwi_fused_b0s.nii.gz -f
chmod 777 topup.sh
./topup.sh > log_topup.txt
# Eddy
LD_LIBRARY_PATH=""
OMP_NUM_THREADS=7
scil_prepare_eddy_command.py ../../split_data/mt_off/dwi.nii.gz ../../split_data/mt_off/dwi.bval ../../split_data/mt_off/dwi.bvec dwi_mean_b0_bet_mask.nii.gz --topup topup_results --out_script --out_prefix dwi_eddy_topup --eddy_cmd eddy_cpu -f
chmod 777 eddy.sh
./eddy.sh > log_eddy.txt
mv dwi_eddy_topup.eddy_rotated_bvecs dwi_eddy_topup.bvec
cp ../../split_data/mt_off/dwi.bval dwi_eddy_topup.bval
cd ../
# MT ON
mkdir mt_on; cd mt_on
# Extract b0 & brain extraction
scil_extract_b0.py ../../split_data/mt_on/dwi.nii.gz ../../split_data/mt_on/dwi.bval ../../split_data/mt_on/dwi.bvec dwi_mean_b0.nii.gz --b0_thr 20 --mean
bet dwi_mean_b0.nii.gz dwi_mean_b0_bet.nii.gz -m -R
# Register rev_b0 to mean_b0
antsRegistrationSyNQuick.sh -d 3 -f dwi_mean_b0.nii.gz -m ../../../raw_data/Anthony_05_19_2022_first_test_WIP_rev_b0_20220519105946_501.nii -o output -t r
mv outputWarped.nii.gz dwi_rev_b0_warped.nii.gz
# Topup
scil_prepare_topup_command.py dwi_mean_b0.nii.gz dwi_rev_b0_warped.nii.gz --out_script --out_prefix topup_results --out_b0s dwi_fused_b0s.nii.gz -f
chmod 777 topup.sh
./topup.sh > log_topup.txt
# Eddy
LD_LIBRARY_PATH=""
OMP_NUM_THREADS=7
scil_prepare_eddy_command.py ../../split_data/mt_on/dwi.nii.gz ../../split_data/mt_on/dwi.bval ../../split_data/mt_on/dwi.bvec dwi_mean_b0_bet_mask.nii.gz --topup topup_results --out_script --out_prefix dwi_eddy_topup --eddy_cmd eddy_cpu -f
chmod 777 eddy.sh
./eddy.sh > log_eddy.txt
mv dwi_eddy_topup.eddy_rotated_bvecs dwi_eddy_topup.bvec
cp ../../split_data/mt_on/dwi.bval dwi_eddy_topup.bval
cd ../../

# N4
mkdir n4; cd n4
# MT OFF
mkdir mt_off; cd mt_off
scil_extract_b0.py ../../eddy_topup/mt_off/dwi_eddy_topup.nii.gz ../../eddy_topup/mt_off/dwi_eddy_topup.bval ../../eddy_topup/mt_off/dwi_eddy_topup.bvec dwi_mean_b0_eddy_topup.nii.gz --mean --b0_thr 20
bet dwi_mean_b0_eddy_topup.nii.gz dwi_mean_b0_eddy_topup_bet.nii.gz -m -R
ln -s -f dwi_mean_b0_eddy_topup_bet_mask.nii.gz mask.nii.gz
N4BiasFieldCorrection -i dwi_mean_b0_eddy_topup.nii.gz -o [dwi_mean_b0_eddy_topup_n4.nii.gz, bias_field_b0.nii.gz] -c [300x150x75x50, 1e-6] > log_N4.txt
scil_apply_bias_field_on_dwi.py ../../eddy_topup/mt_off/dwi_eddy_topup.nii.gz bias_field_b0.nii.gz dwi_eddy_topup_n4.nii.gz --mask mask.nii.gz -f
ln -s -f dwi_eddy_topup_n4.nii.gz dwi.nii.gz
ln -s -f ../../eddy_topup/mt_off/dwi_eddy_topup.bval dwi.bval
ln -s -f ../../eddy_topup/mt_off/dwi_eddy_topup.bvec dwi.bvec
cd ../
# MT ON
mkdir mt_on; cd mt_on
scil_extract_b0.py ../../eddy_topup/mt_on/dwi_eddy_topup.nii.gz ../../eddy_topup/mt_on/dwi_eddy_topup.bval ../../eddy_topup/mt_on/dwi_eddy_topup.bvec dwi_mean_b0_eddy_topup.nii.gz --mean --b0_thr 20
bet dwi_mean_b0_eddy_topup.nii.gz dwi_mean_b0_eddy_topup_bet.nii.gz -m -R
ln -s -f dwi_mean_b0_eddy_topup_bet_mask.nii.gz mask.nii.gz
N4BiasFieldCorrection -i dwi_mean_b0_eddy_topup.nii.gz -o [dwi_mean_b0_eddy_topup_n4.nii.gz, bias_field_b0.nii.gz] -c [300x150x75x50, 1e-6] > log_N4.txt
scil_apply_bias_field_on_dwi.py ../../eddy_topup/mt_on/dwi_eddy_topup.nii.gz bias_field_b0.nii.gz dwi_eddy_topup_n4.nii.gz --mask mask.nii.gz -f
ln -s -f dwi_eddy_topup_n4.nii.gz dwi.nii.gz
ln -s -f ../../eddy_topup/mt_on/dwi_eddy_topup.bval dwi.bval
ln -s -f ../../eddy_topup/mt_on/dwi_eddy_topup.bvec dwi.bvec
cd ../../

# Register MT ON on MT OFF
mkdir registration; cd registration
antsRegistrationSyN.sh -d 3 -n 8 -t r -f ../n4/mt_off/dwi_mean_b0_eddy_topup.nii.gz -m ../n4/mt_on/dwi_mean_b0_eddy_topup.nii.gz -o output
mkdir mt_on
antsApplyTransforms -d 3 -e 3 -i ../n4/mt_on/dwi.nii.gz -r ../n4/mt_off/dwi_mean_b0_eddy_topup.nii.gz -o mt_on/dwi.nii.gz -n NearestNeighbor -t output0GenericAffine.mat outputWarped.nii.gz
cd mt_on/
ln -s -f ../../n4/mt_on/dwi.bval dwi.bval
scil_apply_transform_to_bvecs.py ../../n4/mt_on/dwi.bvec ../output0GenericAffine.mat dwi.bvec -f
scil_extract_b0.py dwi.nii.gz dwi.bval dwi.bvec dwi_mean_b0.nii.gz --mean --b0_thr 20
bet dwi_mean_b0.nii.gz dwi_mean_b0_bet.nii.gz -m -R
ln -s -f dwi_mean_b0_bet_mask.nii.gz mask.nii.gz
cd ../
mkdir mt_off; cd mt_off
ln -s -f ../../n4/mt_off/dwi.nii.gz dwi.nii.gz
ln -s -f ../../n4/mt_off/dwi_mean_b0_eddy_topup_bet_mask.nii.gz mask.nii.gz
ln -s -f ../../n4/mt_off/dwi.bval dwi.bval
ln -s -f ../../n4/mt_off/dwi.bvec dwi.bvec
cd ../../

# DTI
mkdir dti; cd dti
mkdir mt_off; cd mt_off
scil_compute_dti_metrics.py ../../registration/mt_off/dwi.nii.gz ../../registration/mt_off/dwi.bval ../../registration/mt_off/dwi.bvec --mask ../../registration/mt_off/mask.nii.gz -f
cd ../
mkdir mt_on; cd mt_on
scil_compute_dti_metrics.py ../../registration/mt_on/dwi.nii.gz ../../registration/mt_on/dwi.bval ../../registration/mt_on/dwi.bvec --mask ../../registration/mt_on/mask.nii.gz -f
cd ../../

###############################################################
# T1 processing (can be done in || to diffusion processing)
###############################################################
mkdir t1; cd t1
# bet T1
bet ../../raw_data/Anthony_05_19_2022_first_test_WIP_T1W_3D_TFE_0.8mm_20220519105946_301.nii t1_bet.nii.gz -R -m -f 0.3
# denoise T1
scil_run_nlmeans.py t1_bet.nii.gz t1_bet_nlm.nii.gz 1 --mask t1_bet_mask.nii.gz --processes 6 -f
# Segment T1
fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -g -o t1_bet_nlm.nii.gz t1_bet_nlm.nii.gz 

# rename output of T1 segmentation
mv t1_bet_nlm_seg_2.nii.gz wm_mask.nii.gz
mv t1_bet_nlm_seg_1.nii.gz gm_mask.nii.gz
mv t1_bet_nlm_seg_0.nii.gz csf_mask.nii.gz
mv t1_bet_nlm_pve_2.nii.gz wm_map.nii.gz
mv t1_bet_nlm_pve_1.nii.gz gm_map.nii.gz
mv t1_bet_nlm_pve_0.nii.gz csf_map.nii.gz
rm -rf t1_bet_nlm_*

mkdir register_natif; cd register_natif
b0=../../n4/mt_off/dwi_mean_b0_eddy_topup_bet.nii.gz
fa=../../dti/mt_off/fa.nii.gz
t1=../t1_bet_nlm.nii.gz
antsRegistration --dimensionality 3 --float 0\
		 --output [output,outputWarped.nii.gz,outputInverseWarped.nii.gz]\
		 --interpolation Linear --use-histogram-matching 0\
		 --winsorize-image-intensities [0.005,0.995]\
		 --initial-moving-transform [$b0,$t1,1]\
		 --transform Rigid['0.2']\
		 --metric MI[$b0,$t1,1,32,Regular,0.25]\
		 --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
		 --smoothing-sigmas 3x2x1x0\
		 --transform Affine['0.2']\
		 --metric MI[$b0,$t1,1,32,Regular,0.25]\
		 --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
		 --smoothing-sigmas 3x2x1x0\
		 --transform SyN[0.1,3,0]\
		 --metric MI[$b0,$t1,1,32]\
		 --metric CC[$fa,$t1,1,4]\
		 --convergence [50x25x10,1e-6,10] --shrink-factors 4x2x1\
		 --smoothing-sigmas 3x2x1 > log_ants.txt
cd ../

antsApplyTransforms -d 3 -i wm_mask.nii.gz -r ../dti/mt_off/fa.nii.gz -o register_natif/wm_mask.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz
antsApplyTransforms -d 3 -i gm_mask.nii.gz -r ../dti/mt_off/fa.nii.gz -o register_natif/gm_mask.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz
antsApplyTransforms -d 3 -i csf_mask.nii.gz -r ../dti/mt_off/fa.nii.gz -o register_natif/csf_mask.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz
antsApplyTransforms -d 3 -i wm_map.nii.gz -r ../dti/mt_off/fa.nii.gz -o register_natif/wm_map.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz
antsApplyTransforms -d 3 -i gm_map.nii.gz -r ../dti/mt_off/fa.nii.gz -o register_natif/gm_map.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz
antsApplyTransforms -d 3 -i csf_map.nii.gz -r ../dti/mt_off/fa.nii.gz -o register_natif/csf_map.nii.gz -n NearestNeighbor -t register_natif/output0GenericAffine.mat register_natif/output1Warp.nii.gz

cd register_natif
scil_image_math.py convert wm_mask.nii.gz wm_mask.nii.gz --data_type int16 -f
scil_image_math.py convert gm_mask.nii.gz gm_mask.nii.gz --data_type int16 -f
scil_image_math.py convert csf_mask.nii.gz csf_mask.nii.gz --data_type int16 -f
cd ../

# CSD
mkdir csd; cd csd
mkdir mt_on; cd mt_on
scil_compute_ssst_frf.py ../../registration/mt_on/dwi.nii.gz ../../registration/mt_on/dwi.bval ../../registration/mt_on/dwi.bvec frf.txt --mask ../../registration/mt_on/mask.nii.gz --mask_wm ../../t1/register_natif/wm_mask.nii.gz -f
scil_compute_ssst_fodf.py ../../registration/mt_on/dwi.nii.gz ../../registration/mt_on/dwi.bval ../../registration/mt_on/dwi.bvec frf.txt fodf.nii.gz --mask ../../registration/mt_on/mask.nii.gz --sh_order 6 --processes 8 -f
scil_compute_fodf_max_in_ventricles.py fodf.nii.gz ../../dti/mt_on/fa.nii.gz ../../dti/mt_on/md.nii.gz --max_value_output ventricles_fodf_max_value.txt -f
max_value=$(cat ventricles_fodf_max_value.txt)
a_threshold=$(echo 2*${max_value}|bc)
scil_compute_fodf_metrics.py fodf.nii.gz --mask ../../registration/mt_on/mask.nii.gz --at ${a_threshold} --abs_peaks_and_values -f
cd ../
mkdir mt_off; cd mt_off
scil_compute_ssst_frf.py ../../registration/mt_off/dwi.nii.gz ../../registration/mt_off/dwi.bval ../../registration/mt_off/dwi.bvec frf.txt --mask ../../registration/mt_off/mask.nii.gz --mask_wm ../../t1/register_natif/wm_mask.nii.gz -f
scil_compute_ssst_fodf.py ../../registration/mt_off/dwi.nii.gz ../../registration/mt_off/dwi.bval ../../registration/mt_off/dwi.bvec frf.txt fodf.nii.gz --mask ../../registration/mt_off/mask.nii.gz --sh_order 6 --processes 8 -f
scil_compute_fodf_max_in_ventricles.py fodf.nii.gz ../../dti/mt_off/fa.nii.gz ../../dti/mt_off/md.nii.gz --max_value_output ventricles_fodf_max_value.txt -f
max_value=$(cat ventricles_fodf_max_value.txt)
a_threshold=$(echo 2*${max_value}|bc)
scil_compute_fodf_metrics.py fodf.nii.gz --mask ../../registration/mt_off/mask.nii.gz --at ${a_threshold} --abs_peaks_and_values -f
cd ../;

# # Extract shells
# mkdir extracted_shells; cd extracted_shells
# mkdir mt_off; cd mt_off
# scil_extract_dwi_shell.py ../../n4/mt_off/dwi.nii.gz ../../n4/mt_off/dwi.bval ../../n4/mt_off/dwi.bvec 0 dwi_b0.nii.gz dwi_b0.bval dwi_b0.bvec -f
# scil_extract_dwi_shell.py ../../n4/mt_off/dwi.nii.gz ../../n4/mt_off/dwi.bval ../../n4/mt_off/dwi.bvec 1500 dwi_b1500.nii.gz dwi_b1500.bval dwi_b1500.bvec -f
# cd ../
# mkdir mt_on; cd mt_on
# scil_extract_dwi_shell.py ../../n4/mt_on/dwi.nii.gz ../../n4/mt_on/dwi.bval ../../n4/mt_on/dwi.bvec 0 dwi_b0.nii.gz dwi_b0.bval dwi_b0.bvec -f
# scil_extract_dwi_shell.py ../../n4/mt_on/dwi.nii.gz ../../n4/mt_on/dwi.bval ../../n4/mt_on/dwi.bvec 1500 dwi_b1500.nii.gz dwi_b1500.bval dwi_b1500.bvec -f
# cd ../../

# # Powder-average
# mkdir pa; cd pa
# mkdir mt_off; cd mt_off
# python ../../../compute_pa.py ../../extracted_shells/mt_off/dwi_b0.nii.gz dwi_b0_pa.nii.gz
# python ../../../compute_pa.py ../../extracted_shells/mt_off/dwi_b1500.nii.gz dwi_b1500_pa.nii.gz
# cd ../
# mkdir mt_on; cd mt_on
# python ../../../compute_pa.py ../../extracted_shells/mt_on/dwi_b0.nii.gz dwi_b0_pa.nii.gz
# python ../../../compute_pa.py ../../extracted_shells/mt_on/dwi_b1500.nii.gz dwi_b1500_pa.nii.gz
# cd ../..

# mkdir registration; cd registration
# mkdir mt_on
# antsRegistrationSyN.sh -d 3 -n 8 -t r -f ../pa/mt_off/dwi_b0_pa.nii.gz -m ../pa/mt_on/dwi_b0_pa.nii.gz -o output_b0_
# antsRegistrationSyN.sh -d 3 -n 8 -t r -f ../pa/mt_off/dwi_b1500_pa.nii.gz -m ../pa/mt_on/dwi_b1500_pa.nii.gz -o output_b1500_
# mv output_b0_Warped.nii.gz mt_on/dwi_b0_pa.nii.gz
# mv output_b1500_Warped.nii.gz mt_on/dwi_b1500_pa.nii.gz
# mkdir mt_off; cd mt_off
# ln -s -f ../../pa/mt_off/dwi_b0_pa.nii.gz dwi_b0_pa.nii.gz
# ln -s -f ../../pa/mt_off/dwi_b1500_pa.nii.gz dwi_b1500_pa.nii.gz
# ln -s -f ../../n4/mt_off/mask.nii.gz mask.nii.gz
# cd ../../

# MT-contrast
# mkdir mt_contrast; cd mt_contrast
# python ../../compute_mt_contrast.py ../registration/mt_off/dwi_b0_pa.nii.gz ../registration/mt_on/dwi_b0_pa.nii.gz mt_contrast_b0.nii.gz --mask ../registration/mt_off/mask.nii.gz
# python ../../compute_mt_contrast.py ../registration/mt_off/dwi_b1500_pa.nii.gz ../registration/mt_on/dwi_b1500_pa.nii.gz mt_contrast_b1500.nii.gz --mask ../registration/mt_off/mask.nii.gz