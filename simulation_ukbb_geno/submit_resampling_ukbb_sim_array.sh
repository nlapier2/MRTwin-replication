#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-ukbb-sim_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-ukbb-sim.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=16G,highp,h_rt=24:00:00  
#$ -t 1-500:1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

i=$((${SGE_TASK_ID}-1))  # zero-index
bfile=filter4
tld=new_h2c_outdir_sim
temp_dir_name=${tld}/temp/
[ ! -d ${tld} ] && mkdir ${tld}  # make dir if it doesn't exist
[ ! -d ${temp_dir_name} ] && mkdir ${temp_dir_name}  # make dir if it doesn't exist


### SETTINGS ###

num_causal=1000
h2ge=0.2
h2go=0.2
num_trios=1000
num_pcs_confound=20  
num_pcs_regress=10   



### FPR SIMS ###

shuf ${bfile}.fam | head -n 100000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids100k
shuf ${bfile}.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps100k
plink --bfile ${bfile} --extract ${temp_dir_name}${i}.snps100k --keep ${temp_dir_name}${i}.ids100k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink100k
./ProPCA/propca -g ${temp_dir_name}${i}.plink100k -k 20 -o ${temp_dir_name}${i}.pcs100k
#rm ${temp_dir_name}${i}.ids100k ${temp_dir_name}${i}.snps100k 
trio_pct=$(bc -l <<< ${num_trios}'/100000*2')

h2c_arr=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7)
for h2c_ind in {0..7}
do
	this_h2c=${h2c_arr[$h2c_ind]}
	dir_name=${tld}/fpr_h2ce_${this_h2c}_h2co_${this_h2c}/
	[ ! -d ${dir_name} ] && mkdir ${dir_name}  # make dir if it doesn't exist
	echo 'python rand_eff_sim_ukbb.py --bfile '${temp_dir_name}${i}'.plink100k --pcfile '${temp_dir_name}${i}'.pcs100kprojections.txt --num_causal '${num_causal}' --h2ge '${h2ge}' --h2go '${h2go}' --h2eo 0.0 --h2ce '${this_h2c}' --h2co '${this_h2c}' --trio_pct '${trio_pct}' --num_pcs_confound '${num_pcs_confound}' --num_pcs_regress '${num_pcs_regress}' --output '${dir_name}${i}' > '${dir_name}${i}'.sim.log'  > ${dir_name}${i}.sim.log
	python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink100k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo 0.0 --h2ce ${this_h2c} --h2co ${this_h2c} --trio_pct ${trio_pct} --num_pcs_confound ${num_pcs_confound} --num_pcs_regress ${num_pcs_regress} --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
	python post_gwas_snp_select.py --indir ${dir_name} --outdir ${dir_name}post_thresh_5e-8/ --pval_thresh 0.00000005 --num_sims 500
done

#rm ${temp_dir_name}${i}.plink50k* 

#rm ${temp_dir_name}${i}.pcs100k*
#
