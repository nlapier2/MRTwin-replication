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
h2eo=0.05  # for power sims only
num_trios=1000
num_pcs_confound=20  # for FPR sims only
num_pcs_regress=10   # for FPR sims only



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
done


### POWER SIMS ###

trio_pct=$(bc -l <<< ${num_trios}'/100000*2')
dir_name=${tld}/power_ext_100k_h2eo_0.05/
echo 'python rand_eff_sim_ukbb.py --bfile '${temp_dir_name}${i}'.plink100k --pcfile '${temp_dir_name}${i}'.pcs100kprojections.txt --num_causal '${num_causal}' --h2ge '${h2ge}' --h2go '${h2go}' --h2eo '${h2eo}' --h2ce 0.0 --h2co 0.0 --trio_pct '${trio_pct}' --num_pcs_confound 0 --num_pcs_regress 0 --output '${dir_name}${i}' > '${dir_name}${i}'.sim.log'  > ${dir_name}${i}.sim.log
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink100k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound 0 --num_pcs_regress 0 --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
#rm ${temp_dir_name}${i}.plink100k*

shuf ${bfile}.fam | head -n 5000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids5k
shuf ${bfile}.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps5k
plink --bfile ${bfile} --extract ${temp_dir_name}${i}.snps5k --keep ${temp_dir_name}${i}.ids5k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink5k
#rm ${temp_dir_name}${i}.ids5k ${temp_dir_name}${i}.snps5k
trio_pct=$(bc -l <<< ${num_trios}'/5000*2')
dir_name=${tld}/power_ext_5k_h2eo_0.05/
echo 'python rand_eff_sim_ukbb.py --bfile '${temp_dir_name}${i}'.plink5k --pcfile '${temp_dir_name}${i}'.pcs100kprojections.txt --num_causal '${num_causal}' --h2ge '${h2ge}' --h2go '${h2go}' --h2eo '${h2eo}' --h2ce 0.0 --h2co 0.0 --trio_pct '${trio_pct}' --num_pcs_confound 0 --num_pcs_regress 0 --output '${dir_name}${i}' > '${dir_name}${i}'.sim.log'  > ${dir_name}${i}.sim.log
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink5k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound 0 --num_pcs_regress 0 --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
#rm ${temp_dir_name}${i}.plink5k* 

shuf ${bfile}.fam | head -n 10000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids10k
shuf ${bfile}.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps10k
plink --bfile ${bfile} --extract ${temp_dir_name}${i}.snps10k --keep ${temp_dir_name}${i}.ids10k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink10k
#rm ${temp_dir_name}${i}.ids10k ${temp_dir_name}${i}.snps10k
trio_pct=$(bc -l <<< ${num_trios}'/10000*2')
dir_name=${tld}/power_ext_10k_h2eo_0.05/
echo 'python rand_eff_sim_ukbb.py --bfile '${temp_dir_name}${i}'.plink10k --pcfile '${temp_dir_name}${i}'.pcs100kprojections.txt --num_causal '${num_causal}' --h2ge '${h2ge}' --h2go '${h2go}' --h2eo '${h2eo}' --h2ce 0.0 --h2co 0.0 --trio_pct '${trio_pct}' --num_pcs_confound 0 --num_pcs_regress 0 --output '${dir_name}${i}' > '${dir_name}${i}'.sim.log'  > ${dir_name}${i}.sim.log
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink10k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound 0 --num_pcs_regress 0 --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
#rm ${temp_dir_name}${i}.plink10k* 

shuf ${bfile}.fam | head -n 20000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids20k
shuf ${bfile}.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps20k
plink --bfile ${bfile} --extract ${temp_dir_name}${i}.snps20k --keep ${temp_dir_name}${i}.ids20k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink20k
#rm ${temp_dir_name}${i}.ids20k ${temp_dir_name}${i}.snps20k
trio_pct=$(bc -l <<< ${num_trios}'/20000*2')
dir_name=${tld}/power_ext_20k_h2eo_0.05/
echo 'python rand_eff_sim_ukbb.py --bfile '${temp_dir_name}${i}'.plink20k --pcfile '${temp_dir_name}${i}'.pcs100kprojections.txt --num_causal '${num_causal}' --h2ge '${h2ge}' --h2go '${h2go}' --h2eo '${h2eo}' --h2ce 0.0 --h2co 0.0 --trio_pct '${trio_pct}' --num_pcs_confound 0 --num_pcs_regress 0 --output '${dir_name}${i}' > '${dir_name}${i}'.sim.log'  > ${dir_name}${i}.sim.log
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink20k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound 0 --num_pcs_regress 0 --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
#rm ${temp_dir_name}${i}.plink20k* 

shuf ${bfile}.fam | head -n 50000 | awk '{print $1,$2}' > ${temp_dir_name}${i}.ids50k
shuf ${bfile}.bim | head -n 22000 | awk '{print $2}' > ${temp_dir_name}${i}.snps50k
plink --bfile ${bfile} --extract ${temp_dir_name}${i}.snps50k --keep ${temp_dir_name}${i}.ids50k --geno 0.0 --make-bed --out ${temp_dir_name}${i}.plink50k
#rm ${temp_dir_name}${i}.ids50k ${temp_dir_name}${i}.snps50k
trio_pct=$(bc -l <<< ${num_trios}'/50000*2')
dir_name=${tld}/power_ext_50k_h2eo_0.05/
echo 'python rand_eff_sim_ukbb.py --bfile '${temp_dir_name}${i}'.plink50k --pcfile '${temp_dir_name}${i}'.pcs100kprojections.txt --num_causal '${num_causal}' --h2ge '${h2ge}' --h2go '${h2go}' --h2eo '${h2eo}' --h2ce 0.0 --h2co 0.0 --trio_pct '${trio_pct}' --num_pcs_confound 0 --num_pcs_regress 0 --output '${dir_name}${i}' > '${dir_name}${i}'.sim.log'  > ${dir_name}${i}.sim.log
python rand_eff_sim_ukbb.py --bfile ${temp_dir_name}${i}.plink50k --pcfile ${temp_dir_name}${i}.pcs100kprojections.txt --num_causal ${num_causal} --h2ge ${h2ge} --h2go ${h2go} --h2eo ${h2eo} --h2ce 0.0 --h2co 0.0 --trio_pct ${trio_pct} --num_pcs_confound 0 --num_pcs_regress 0 --output ${dir_name}${i} >> ${dir_name}${i}.sim.log
#rm ${temp_dir_name}${i}.plink50k* 

#rm ${temp_dir_name}${i}.pcs100k*
#
