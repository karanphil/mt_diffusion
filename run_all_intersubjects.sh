#! /bin/bash
cd /home/pkaran/Samsung/data/MT_Diffusion/myelo_inferno/output;
all_dirs=$(ls -d sub*);
specs="1_degrees_bins_0.5_FA_thr_NuFo_False";
subjects=("");
last_sub="None";
for dir in $all_dirs;
    do IFS=',';
    IFS='_';
    read -a strarr <<< $dir;
    IFS=',';
    subject=${strarr[0]};
    # sub_dirs=$(ls -d ${subject}*);
    if [ $last_sub != $subject ];
        then subjects+=,$subject;
        last_sub=$subject
    fi;
done;
python ~/source/mt_diffusion/plot_inter_subject.py intersubject/all_subjects_${specs}.png${subjects};