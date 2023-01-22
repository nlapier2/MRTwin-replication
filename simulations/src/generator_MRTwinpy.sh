#!/bin/sh

source /u/local/Modules/default/init/modules.sh
module load R/4.0.2
module load gcc
module load anaconda3

oFile=../data/ # change here

repeat=1
snpfilter="GWAS"


codefile=generator.R
odir="results_MR"
tflag=1 # tflag = 1 means perform GWAS filter using external data, otherwise all the SNPs will be used for the analysis
SGE_TASK_ID=1 # SGE_TASK_ID correspond to the row of the parameters set (1-based) in the demo.param.txt file
permute=100 # number of permutation used for MR-Twin

realD=$SGE_TASK_ID
realD=$((realD-1))
echo "realD is $realD"
RowIndex=$((realD/repeat))
RowIndex=$((RowIndex+2))
Seed=$((realD%repeat))
Seed=$((Seed+1))
echo "Seed is ${Seed}, RowIndex is ${RowIndex}"
paramPath='../param/demo.param.txt'
oPath="../${odir}_${snpfilter}/"
if [ ! -d $oPath ]; then
  mkdir $oPath
fi
echo "This experiments run power and calibration testing on MR-duo, sib" | tee ${oPath}readme.setting.log
code=${codefile}
# mapfile -t < /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/MR/sriram_sim/results_weakIV_fine_grid/polygenic_brumpton_sim_fine_grid.txt
var=$(awk "NR==${RowIndex}" ${paramPath})
IFS=' ' read nsim m mcausal n h2 se so gamma_ue gamma_uo ce ext_n ext_model trio_model efst fst <<< $var
nsim=3
echo "${code} ${Seed} $nsim $m $mcausal $n $h2 $se $so $gamma_ue $gamma_uo $ce $ext_n $ext_model $trio_model $efst $fst $snpfilter $realD $tflag" 
Rscript ${code} ${Seed} $nsim $m $mcausal $n $h2 $se $so $gamma_ue $gamma_uo $ce $ext_n $ext_model $trio_model $efst $fst  $snpfilter $realD $tflag $oFile

echo "nsim: ${nsim}, m: $m, mcausal: $mcausal, n: $n, h2: $h2, se: $se, so: $so,  ue: $gamma_ue, uo: $gamma_uo, ce: $ce, ext_n: $ext_n, ext_model: $ext_model, trio_model: $trio_model, efst: $efst, fst: $fst, snpfilter:  $snpfilter" | tee ${oPath}${RowIndex}_${Seed}_$((realD)).MR.log

dataPath="$oFile/"
for (( rep=1; rep<=$nsim; rep++ )); do
  datafile="${realD}.${Seed}.${rep}"
  # scripts=/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/MR/MR-Twin/demo/scripts/
  scripts=../MR-Twin/
  echo "*** rep ${rep} ***" 

  # python ${scripts}/run_mrtwin_real_data.py --include_duo --include_sib --num_twins 100 --basename ${dataPath}/${datafile} >> ${oPath}${RowIndex}_${Seed}_$((realD)).MR.log
  python ${scripts}/run_mrtwin_real_data.py --num_twins $permute --basename ${dataPath}/${datafile} | tee -a ${oPath}${RowIndex}_${Seed}_$((realD)).MR.log

  echo "removed the tmp files from ${dataPath}/${datafile}"

  rm ${dataPath}/${datafile}.betahatEO
  rm ${dataPath}/${datafile}.pa1geno
  rm ${dataPath}/${datafile}.pa2geno
  rm ${dataPath}/${datafile}.childgeno
  rm ${dataPath}/${datafile}.childpheno
  rm ${dataPath}/${datafile}.child2geno
  rm ${dataPath}/${datafile}.child2pheno

  echo "finish removal"
done
# Rscript ${code} ${Seed} $nsim $m $mcausal $n $h2 $se $so $gamma_ue $gamma_uo $ce $ext_n $ext_model $trio_model $efst $fst $snpfilter > ${oPath}${RowIndex}_${Seed}_$((realD))_M100.log
# fi
