#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-ukbb-sim_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-ukbb-sim.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=16G,highp,h_rt=24:00:00  
#$ -t 1-200:1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

i=$((${SGE_TASK_ID}-1))  # zero-index
#tld=outdir_sim
tld=sample_size_fpr_scaling
temp_dir_name=${tld}/temp/
[ ! -d ${tld} ] && mkdir ${tld}  # make dir if it doesn't exist
[ ! -d ${temp_dir_name} ] && mkdir ${temp_dir_name}  # make dir if it doesn't exist


### SETTINGS ###

num_causal=1000
h2ge=0.2
h2go=0.2
this_h2c=0.3
h2eo=0.05  # for power sims only
num_trios=1000
num_pcs_confound=20  # for FPR sims only
num_pcs_regress=10   # for FPR sims only


# 100k

shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.fam | head -n 100000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids100k
shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps100k
plink --bfile /u/project/sriram/ukbiobank/data/geno/cal/filter4 --extract ${temp_dir_name}${i}.snps100k --keep ${temp_dir_name}${i}.ids100k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink100k
./ProPCA/propca -g ${temp_dir_name}${i}.plink100k -k 20 -o ${temp_dir_name}${i}.pcs100k
trio_pct=$(bc -l <<< ${num_trios}'/100000*2')

dir_name=${tld}/fpr_100k/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink100k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo 0.0 --h2ce ${this_h2c} --h2co ${this_h2c} --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000

dir_name=${tld}/power_ext_100k_h2eo_0.05/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink100k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000
rm ${temp_dir_name}${i}.plink100k*
rm ${temp_dir_name}${i}.pcs100k*
rm ${dir_name}${i}*plink*



# 5k

shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.fam | head -n 5000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids5k
shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps5k
plink --bfile /u/project/sriram/ukbiobank/data/geno/cal/filter4 --extract ${temp_dir_name}${i}.snps5k --keep ${temp_dir_name}${i}.ids5k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink5k
./ProPCA/propca -g ${temp_dir_name}${i}.plink5k -k 20 -o ${temp_dir_name}${i}.pcs5k
trio_pct=$(bc -l <<< ${num_trios}'/5000*2')

dir_name=${tld}/fpr_5k/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink5k --pcfile ${temp_dir_name}${i}.pcs5kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo 0.0 --h2ce ${this_h2c} --h2co ${this_h2c} --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000

dir_name=${tld}/power_ext_5k_h2eo_0.05/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink5k --pcfile ${temp_dir_name}${i}.pcs5kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000
rm ${temp_dir_name}${i}.plink5k* 
rm ${temp_dir_name}${i}.pcs5k*
rm ${dir_name}${i}*plink*



# 10k

shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.fam | head -n 10000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids10k
shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps10k
plink --bfile /u/project/sriram/ukbiobank/data/geno/cal/filter4 --extract ${temp_dir_name}${i}.snps10k --keep ${temp_dir_name}${i}.ids10k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink10k
./ProPCA/propca -g ${temp_dir_name}${i}.plink10k -k 20 -o ${temp_dir_name}${i}.pcs10k
trio_pct=$(bc -l <<< ${num_trios}'/10000*2')

dir_name=${tld}/fpr_10k/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink10k --pcfile ${temp_dir_name}${i}.pcs10kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo 0.0 --h2ce ${this_h2c} --h2co ${this_h2c} --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000

dir_name=${tld}/power_ext_10k_h2eo_0.05/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink10k --pcfile ${temp_dir_name}${i}.pcs10kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000
rm ${temp_dir_name}${i}.plink10k* 
rm ${temp_dir_name}${i}.pcs10k*
rm ${dir_name}${i}*plink*



# 20k

shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.fam | head -n 20000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids20k
shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps20k
plink --bfile /u/project/sriram/ukbiobank/data/geno/cal/filter4 --extract ${temp_dir_name}${i}.snps20k --keep ${temp_dir_name}${i}.ids20k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink20k
./ProPCA/propca -g ${temp_dir_name}${i}.plink20k -k 20 -o ${temp_dir_name}${i}.pcs20k
trio_pct=$(bc -l <<< ${num_trios}'/20000*2')

dir_name=${tld}/fpr_20k/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink20k --pcfile ${temp_dir_name}${i}.pcs20kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo 0.0 --h2ce ${this_h2c} --h2co ${this_h2c} --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000

dir_name=${tld}/power_ext_20k_h2eo_0.05/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink20k --pcfile ${temp_dir_name}${i}.pcs20kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000
rm ${temp_dir_name}${i}.plink20k* 
rm ${temp_dir_name}${i}.pcs20k*
rm ${dir_name}${i}*plink*



# 50k

shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.fam | head -n 50000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids50k
shuf /u/project/sriram/ukbiobank/data/geno/cal/filter4.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps50k
plink --bfile /u/project/sriram/ukbiobank/data/geno/cal/filter4 --extract ${temp_dir_name}${i}.snps50k --keep ${temp_dir_name}${i}.ids50k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink50k
./ProPCA/propca -g ${temp_dir_name}${i}.plink50k -k 20 -o ${temp_dir_name}${i}.pcs50k
trio_pct=$(bc -l <<< ${num_trios}'/50000*2')

dir_name=${tld}/fpr_50k/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink50k --pcfile ${temp_dir_name}${i}.pcs50kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo 0.0 --h2ce ${this_h2c} --h2co ${this_h2c} --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000

dir_name=${tld}/power_ext_50k_h2eo_0.05/
[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink50k --pcfile ${temp_dir_name}${i}.pcs50kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 1000
rm ${temp_dir_name}${i}.plink50k* 
rm ${temp_dir_name}${i}.pcs50k*
rm ${dir_name}${i}*plink*
#
