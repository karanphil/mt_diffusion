#! /bin/bash
cd ~/Samsung/data/MT_Diffusion/myelo_inferno/DTI_metrics;
all_dirs=$(ls);

output_choice=$1;

for dir in $all_dirs;
    do
    if [[ $output_choice == "ratios" ]];
        then cd ../output_ratios;
        mt_path="../ihMT/${dir}/${dir}__MTR_warped.nii.gz";
        ihmt_path="../ihMT/${dir}/${dir}__ihMTR_warped.nii.gz";
    elif [[ $output_choice == "sats" ]];
        then cd ../output_sats;
        mt_path="../ihMT/${dir}/${dir}__MTsat_warped.nii.gz";
        ihmt_path="../ihMT/${dir}/${dir}__ihMTsat_warped.nii.gz";
    else
        echo "Invalid choice.";
        exit 0;
    fi;
    echo $dir;
    python ~/Research/source/mt_diffusion/average_mt_ihmt_by_angle_dti.py ../DTI_metrics/${dir}/${dir}__dti_evecs_v1.nii.gz ../DTI_metrics/${dir}/${dir}__dti_fa.nii.gz $mt_path $ihmt_path ../wm_mask/${dir}/${dir}__wm_mask.nii.gz ../FODF_metrics/${dir}/nufo.nii.gz $output_choice;
done;
