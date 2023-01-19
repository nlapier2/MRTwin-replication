#!/bin/bash
#$ -cwd
#$ -l h_data=4G,h_rt=0:30:00,highp
#$ -t 1-12:1
#$ -N job_prune
#$ -o stdout_prune

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

exposure_list=prune_exposures.txt
mapfile EXPOS -d $'\n' -t < $exposure_list
exposure=${EXPOS[$((SGE_TASK_ID-1))]}
exposure=${exposure//$' '}
exposure=${exposure//$'\n'}

python scripts/ld_prune.py --gwas gwas/gwas_combine_standard/${exposure}.gwas.pheno.glm.linear --bfile /u/project/sgss/UKBB/data/cal/filter4 --gwas_thresh 0.00000005 --prune_thresh 0.1 --temp_basename ${exposure} --out pruned/${exposure}.prune
python scripts/generate_cor_mat.py --bfile /u/project/sgss/UKBB/data/cal/filter4 --snps pruned/${exposure}.prune.in --out pruned/${exposure}.prune
#
