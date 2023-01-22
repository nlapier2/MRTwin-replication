#!/bin/sh

source /u/local/Modules/default/init/modules.sh
module load R/4.0.2
module load gcc

SGE_TASK_ID=1 # SGE_TASK_ID correspond to the row of the parameters set (1-based) in the demo.param.txt file

repeat=1
snpfilter="GWAS"

codefile=benchmark.R
odir="results_benchmark"

realD=$SGE_TASK_ID
echo $realD
realD=$((realD-1))
# echo "realD is $realD"
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
echo "This experiment serves as benchmarking different methods" | tee ${oPath}readme.setting.log
code=${codefile}
# mapfile -t < /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/MR/sriram_sim/results_weakIV_fine_grid/polygenic_brumpton_sim_fine_grid.txt
var=$(awk "NR==${RowIndex}" ${paramPath})
IFS=' ' read nsim m mcausal n h2 se so gamma_ue gamma_uo ce ext_n ext_model trio_model efst fst <<< $var
# echo "${MAPFILE[$((SGE_TASK_ID-1))]}"
echo "${code} ${Seed} $nsim $m $mcausal $n $h2 $se $so $gamma_ue $gamma_uo $ce $ext_n $ext_model $trio_model $efst $fst $snpfilter"
# Rscript ${code} ${Seed} $nsim $m $mcausal $n $h2 $se $so $gamma_ue $gamma_uo $ce $ext_n $ext_model $trio_model $efst $fst  $snpfilter
Rscript ${code} ${Seed} $nsim $m $mcausal $n $h2 $se $so $gamma_ue $gamma_uo $ce $ext_n $ext_model $trio_model $efst $fst $snpfilter | tee -a ${oPath}${RowIndex}_${Seed}_$((realD))_M100.log
# fi
