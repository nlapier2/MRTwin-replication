# Instructions

1. We are not able to include the PLINK bfile that we used in our simulations because we used genotypes from the UK Biobank project. Please follow the instructions outlined in the paper to reproduce this bfile. This essentially involves isolating all unrelated White British individuals from the UK Biobank. 

2. Replace the bfile name (which is "filter4" in ours) in "submit\_resampling\_ukbb\_sim\_array.sh" with the one you created.

3. Run "submit\_resampling\_ukbb\_sim\_array.sh", possibly modifying it if needed to work with your job scheduler. It is not recommended to run this on a laptop, as it will consume a great deal of time and storage (several hundred GB).

4. Run "submit\_eval\_main.sh" and examine the results.

5. Optionally, if you would like, you can run duo and sib modes by using the --include\_duo and/or --include\_sib flags in the calls to eval\_sim\_mrtrio.py in the submission script for step 3. Please be advised that these take much longer to run, especially sib mode.

