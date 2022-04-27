#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-eval_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-eval.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=8G,highp,h_rt=16:00:00  
#$ -t 1-13:1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

dir_arr=(new_h2c_outdir_sim/fpr_h2ce_0.0_h2co_0.0/ new_h2c_outdir_sim/fpr_h2ce_0.1_h2co_0.1/ new_h2c_outdir_sim/fpr_h2ce_0.2_h2co_0.2/ new_h2c_outdir_sim/fpr_h2ce_0.3_h2co_0.3/ new_h2c_outdir_sim/fpr_h2ce_0.4_h2co_0.4/ new_h2c_outdir_sim/fpr_h2ce_0.5_h2co_0.5/ new_h2c_outdir_sim/fpr_h2ce_0.6_h2co_0.6/ new_h2c_outdir_sim/fpr_h2ce_0.7_h2co_0.7/ new_h2c_outdir_sim/power_ext_5k_h2eo_0.05/ new_h2c_outdir_sim/power_ext_10k_h2eo_0.05/ new_h2c_outdir_sim/power_ext_20k_h2eo_0.05/ new_h2c_outdir_sim/power_ext_50k_h2eo_0.05/ new_h2c_outdir_sim/power_ext_100k_h2eo_0.05/) 

index=$((${SGE_TASK_ID}-1))
dir_name=${dir_arr[$index]}
num_sims=500
num_twins=100

python eval_sim_mrtrio.py --dir_name ${dir_name} --method ivw --num_sims ${num_sims} --num_twins ${num_twins} > ${dir_name}results_ivw_mrtrio_corr.txt
Rscript brumpton.R ${dir_name} ${num_sims} > ${dir_name}results_brumpton.txt
python eval_sim_mrtrio.py --dir_name ${dir_name} --method egger --num_sims ${num_sims} --num_twins ${num_twins} > ${dir_name}results_egger_mrtrio.txt
python eval_sim_mrtrio.py --dir_name ${dir_name} --method median --num_sims ${num_sims} --num_twins ${num_twins} > ${dir_name}results_median_mrtrio.txt
python eval_sim_mrtrio.py --dir_name ${dir_name} --method mode --num_sims ${num_sims} --num_twins ${num_twins} > ${dir_name}results_mode_mrtrio.txt
#
