#!/bin/sh

mapfile -t < params.txt
source /u/local/Modules/default/init/modules.sh
SGE_TASK_ID=1
module load R
module load gcc
Rscript mrtrio_controlPC.R "${MAPFILE[$((SGE_TASK_ID))]}" 
# fi
