#!/bin/bash
#$ -cwd
#$ -l h_data=16G,h_rt=0:30:00,highp
#$ -t 1-144:1
#$ -N job_trait_pairs1
#$ -o stdout_trait_pairs_1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

intermed_dir=zzz_intermediate_ext/
results_dir=zzz_results_ext/
trio_dir=/u/project/sriram/boyang19/trioExtract_orig/
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

python scripts/make_mrtwin_infiles.py --trio_file ${trio_dir}trio_conserve_white_British.txt --path_to_pheno ${trio_dir}traitlist/combine_standard/ --ext_bfile /u/project/sgss/UKBB/data/cal/filter4 --child_bfile ${trio_dir}trios/British_child --p1_bfile ${trio_dir}trios/P1_WB --p2_bfile ${trio_dir}trios/P2_WB --exposure_assoc gwas/gwas_combine_standard/${exposure}.gwas.pheno.glm.linear --outcome_assoc gwas/gwas_combine_standard/${outcome}.gwas.pheno.glm.linear --exposure_name ${exposure} --outcome_name ${outcome} --outname ${intermed_dir}${pairname} --snps pruned/${exposure}.prune.in

outname=${results_dir}/${pairname}_results.txt
python scripts/run_mrtwin_real_data.py --basename ${intermed_dir}${pairname} --num_twins 1000 > ${outname}
unset R_HOME
Rscript scripts/run_mr.R ${intermed_dir}${pairname}.betahatEO pruned/${exposure}.prune.ld >> ${outname}
Rscript scripts/brumpton_filter.R ${intermed_dir}${pairname} ${num_sims} >> ${outname}
#
