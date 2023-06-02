#! /bin/bash
output_choice=$1;
all_dirs=(sub*);
if [[ $output_choice = "ratios" ]];
    then rs="r";
elif [[ $output_choice = "sats" ]];
    then rs="sat";
else
    echo "Invalid choice.";
    exit 0;
fi;
specs="WB_single_fiber_mt${rs}_ihmt${rs}_results_0.5_fa_thr_1_bin_width";
specs_cr="WB_single_fiber_corrected_mt${rs}_ihmt${rs}_results_0.5_fa_thr_1_bin_width";
subjects=("");
last_sub="None";
default_IFS=$IFS
for dir in "${all_dirs[@]}";
    do IFS='_';
    read -a strarr <<< $dir;
    IFS=$default_IFS;
    subject=${strarr[0]};
    if [ $last_sub != $subject ];
        then subjects+=" "$subject;
        last_sub=$subject
    fi;
done;
python ~/source/mt_diffusion/scripts/plot_intersubject_corrected.py intersubject/all_subjects_${specs}.png --in_name results_txt/${specs}.txt --input_type ${output_choice} --in_subjects ${subjects};