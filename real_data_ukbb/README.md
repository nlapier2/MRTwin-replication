Instructions:

1. We are not able to include the PLINK bfile that we used in our simulations because we used genotypes from the UK Biobank project. Please follow the instructions outlined in the paper to reproduce this bfile. This essentially involves isolating all unrelated White British individuals from the UK Biobank. Similarly, the phenotype files and trio identitites cannot be publicly released and must be re-generated according to the instructions in the paper. In the scripts below, the bfile is called "filter4", the directory with the standardized phenotypes is "pheno\_combine\_standard/", and the trio files are located in "trioExtract\_orig/". You will have to replace the names in the submission scripts with the files/directories you generate, or name them the same thing.

2. Run submit\_gwas.sh to get the betas. Store in some directory -- in this example, gwas/gwas\_combine\_standard/. 

3. Then run submit\_ld\_prune.sh to extract and prune significant SNPs. Store in some directory -- in this example, pruned/.

4. Then run submit\_mr\_ext.sh MR submit scripts.


