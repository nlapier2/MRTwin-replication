#!/bin/sh
#$ -cwd
#$ -l h_data=2G,h_rt=0:10:00,highp
#$ -t 1-2:1
#$ -N job-runGWAS
#$ -o stdout-runGWAS

source ~/.bash_profile
source /u/local/Modules/default/init/modules.sh
path=PATH ##### CHANGE

#SETTINGS
ps_1_null=$path/ps_1_null
ps_2_null=$path/ps_2_null
ps_2_alter=$path/ps_2_alter

setting=SETTING ##### CHANGE: Choose one from above i.e. $ps_1_null 

bfile=$setting/external/external.X
outdir=gwas/gwas_combine_standard/
dircp=$setting/external/
phenofile=pheno.txt
mapfile -t < $phenofile
pheno=${MAPFILE[$((SGE_TASK_ID-1))]}
#plink2 --bfile ${bfile} --linear hide-covar --variance-standardize --pheno  ${dircp}${pheno}.pheno --out ${outdir}${pheno}.gwas --covar ${dircp}${pheno}.covar
plink2 --bfile ${bfile} --linear allow-no-covars --variance-standardize --pheno  ${dircp}${pheno}.pheno --out ${outdir}${pheno}.gwas
#
