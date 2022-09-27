#!/bin/bash
#$ -cwd
#$ -l h_data=2G,h_rt=0:30:00,highp
#$ -t 1-1:1
#$ -N job_trait_pairs1
#$ -o stdout_trait_pairs_1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile


path=PATH ##### CHANGE

#SETTINGS
ps_1_null=$path/ps_1_null
ps_2_null=$path/ps_2_null
ps_2_alter=$path/ps_2_alter

setting=SETTING ##### CHANGE: Choose one from above i.e. $ps_1_null


intermed_dir=intermediate_ext/
results_dir=results_ext
trio_dir=$setting/ps_2_alter/
exposure_list=full_exposures.txt
outcome_list=full_outcomes.txt

mapfile EXPOS -d $'\n' -t < $exposure_list
mapfile OUTCOME_V -d $'\n' -t < $outcome_list
exposure=${EXPOS[$((SGE_TASK_ID-1))]}
outcome=${OUTCOME_V[$((SGE_TASK_ID-1))]}
exposure=${exposure//$' '}
exposure=${exposure//$'\n'}
outcome=${outcome//$'\n'}
outcome=${outcome//$' '}
pairname=${exposure}---${outcome}

python3 scripts/make_mrtwin_infiles.py --trio_file ${trio_dir}trios/trio.txt --path_to_pheno ${trio_dir}trios/ --ext_bfile ${trio_dir}external/external.X --child_bfile ${trio_dir}trios/trio.child --p1_bfile ${trio_dir}trios/trio.father --p2_bfile ${trio_dir}trios/trio.mother --exposure_assoc gwas/gwas_combine_standard/${exposure}.gwas.pheno.glm.linear --outcome_assoc gwas/gwas_combine_standard/${outcome}.gwas.pheno.glm.linear --exposure_name ${exposure} --outcome_name ${outcome} --outname ${intermed_dir}${pairname} --snps pruned/${exposure}.prune.in
outname=${results_dir}/${pairname}_results.txt
cat "ivw results:" > ${outname}
python scripts/run_mrtwin_real_data.py --basename ${intermed_dir}${pairname} > ${outname}
