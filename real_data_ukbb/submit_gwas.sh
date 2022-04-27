#!/bin/sh
#$ -cwd
#$ -l h_data=16G,h_rt=3:30:00,highp
#$ -t 1-12:1
#$ -N job-runGWAS
#$ -o stdout-runGWAS

bfile=filter4
outdir=gwas_combine_standard/
dircp=pheno_combine_standard/
phenofile=pheno.txt
mapfile -t < $phenofile
pheno=${MAPFILE[$((SGE_TASK_ID-1))]}
plink2 --bfile ${bfile} --linear hide-covar --variance-standardize --pheno  ${dircp}${pheno}.pheno --out ${outdir}${pheno}.gwas
#
