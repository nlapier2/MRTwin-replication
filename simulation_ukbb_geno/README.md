# Replication instructions for the UK Biobank Simulations

This directory contains scripts and instructions for replicating the UK Biobank simulation results in the MR-Twin main text.


### Instructions for the main FPR figure

1. We are not able to include the PLINK bfile that we used in our simulations because we used genotypes from the UK Biobank project. Please follow the instructions outlined in the paper to reproduce this bfile. This essentially involves isolating all unrelated White British individuals from the UK Biobank. 

2. Replace the bfile name (which is "filter4" in ours) in "submit\_resampling\_ukbb\_sim\_array.sh" with the one you created.

3. Run "submit\_resampling\_ukbb\_sim\_array.sh", possibly modifying it if needed to work with your job scheduler. It is not recommended to run this on a laptop, as it will consume a great deal of time and storage (several hundred GB).

4. Run "submit\_eval\_main.sh" and examine the results.


### Instructions for other figures

For replicating the other figures, you can follow the steps above but swap out the scripts in the last two steps with the ones below.

Replicating duo and sib mode results as well as results using different MR methods with MR-Twin: see "submit\_eval\_duo\_sib\_stats.sh". In practice, duo and sib modes take much longer to run, so it is probably preferable to parallelize jobs and then combine the results into a single file later.

Replicating the MR-Twin power and FPR results under varying sample sizes: see "submit\_sample\_size\_fpr\_scaling.sh" and "submit\_eval\_sample\_size\_fpr.sh".

