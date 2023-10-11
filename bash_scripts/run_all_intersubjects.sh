#! /bin/bash
output_choice=$1;
bin_width=$2;
if [[ $output_choice = "ratios" ]];
    then cd ~/Samsung/data/MT_Diffusion/myelo_inferno/output_ratios;
elif [[ $output_choice = "sats" ]];
    then cd ~/Samsung/data/MT_Diffusion/myelo_inferno/output_sats;
else
    echo "Invalid choice.";
    exit 0;
fi;
all_dirs=$(ls -d sub*);
specs="${bin_width}_degrees_bins_0.5_FA_thr_NuFo_False";
subjects=("");
last_sub="None";
default_IFS=$IFS
for dir in $all_dirs;
    do IFS='_';
    read -a strarr <<< $dir;
    IFS=$default_IFS;
    subject=${strarr[0]};
    # sub_dirs=$(ls -d ${subject}*);
    if [ $last_sub != $subject ];
        then subjects+=" "$subject;
        last_sub=$subject
    fi;
done;
python ~/Research/source/mt_diffusion/plot_inter_subject.py intersubject/all_subjects_${specs}.png $output_choice $bin_width ${subjects};