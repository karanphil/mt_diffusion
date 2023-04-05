#! /bin/bash
# cd /home/pkaran/Samsung/data/MT_Diffusion/myelo_inferno/output_ratios;
cd /home/pkaran/Samsung/data/MT_Diffusion/myelo_inferno/output_sats;
all_dirs=$(ls -d sub*);
specs="1_degrees_bins_0.5_FA_thr_NuFo_False";
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
python ~/source/mt_diffusion/plot_inter_subject.py intersubject/all_subjects_${specs}.png ${subjects};