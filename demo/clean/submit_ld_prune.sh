#!/bin/bash
#$ -cwd
#$ -l h_data=4G,h_rt=0:30:00,highp
#$ -t 1-12:1
#$ -N job_prune
#$ -o stdout_prune

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

path=PATH ##### CHANGE

#SETTINGS
ps_1_null=$path/ps_1_null
ps_2_null=$path/ps_2_null
ps_2_alter=$path/ps_2_alter

setting=SETTING ##### CHANGE: Choose one from above i.e. $ps_1_null


exposure_list=prune_exposures.txt
mapfile EXPOS -d $'\n' -t < $exposure_list
exposure=${EXPOS[$((SGE_TASK_ID-1))]}
exposure=${exposure//$' '}
exposure=${exposure//$'\n'}

python3 scripts/ld_prune.py --gwas gwas/gwas_combine_standard/${exposure}.gwas.pheno.glm.linear --bfile $setting/external/external.X --gwas_thresh 0.01 --prune_thresh 0.1 --temp_basename ${exposure} --out pruned/${exposure}.prune
