#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-eval_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-eval.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=8G,highp,h_rt=16:00:00  
#$ -t 1-20:1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

dir_arr=(sample_size_fpr_scaling/fpr_100k/ sample_size_fpr_scaling/fpr_5k/ sample_size_fpr_scaling/fpr_10k/ sample_size_fpr_scaling/fpr_20k/ sample_size_fpr_scaling/fpr_50k/ sample_size_fpr_scaling/power_ext_5k_h2eo_0.05/ sample_size_fpr_scaling/power_ext_10k_h2eo_0.05/ sample_size_fpr_scaling/power_ext_20k_h2eo_0.05/ sample_size_fpr_scaling/power_ext_50k_h2eo_0.05/ sample_size_fpr_scaling/power_ext_100k_h2eo_0.05/ sample_size_fpr_scaling/fpr_100k/post_thresh_5e-8/ sample_size_fpr_scaling/fpr_5k/post_thresh_5e-8/ sample_size_fpr_scaling/fpr_10k/post_thresh_5e-8/ sample_size_fpr_scaling/fpr_20k/post_thresh_5e-8/ sample_size_fpr_scaling/fpr_50k/post_thresh_5e-8/ sample_size_fpr_scaling/power_ext_5k_h2eo_0.05/post_thresh_5e-8/ sample_size_fpr_scaling/power_ext_10k_h2eo_0.05/post_thresh_5e-8/ sample_size_fpr_scaling/power_ext_20k_h2eo_0.05/post_thresh_5e-8/ sample_size_fpr_scaling/power_ext_50k_h2eo_0.05/post_thresh_5e-8/ sample_size_fpr_scaling/power_ext_100k_h2eo_0.05/post_thresh_5e-8/)

index=$((${SGE_TASK_ID}-1))
dir_name=${dir_arr[$index]}
num_sims=200
num_twins=100

python eval_sim_mrtrio.py --dir_name ${dir_name} --method ivw --num_sims ${num_sims} --num_twins ${num_twins} > ${dir_name}results_ivw_mrtrio.txt
#
