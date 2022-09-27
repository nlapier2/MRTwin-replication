path=full_example/files

#SETTINGS
ps_1_null=$path/ps_1_null
ps_2_null=$path/ps_2_null
ps_2_alter=$path/ps_2_alter

setting=$ps_2_null  ##### CHANGE: Choose one from above i.e. $ps_1_null 

bfile=$setting/external/external.X
trio_dir=$setting/
outdir=full_example/
exposure=exposure
outcome=outcome
pairname=${exposure}---${outcome}
phenofile1=${setting}/external/exposure.pheno
phenofile2=${setting}/external/outcome.pheno

plink2 --bfile ${bfile} --linear --variance-standardize --pheno  ${phenofile1} --out ${outdir}${exposure}.gwas
plink2 --bfile ${bfile} --linear --variance-standardize --pheno  ${phenofile2} --out ${outdir}${outcome}.gwas

python3 scripts/ld_prune.py --gwas ${outdir}${exposure}.gwas.PHENO1.glm.linear --bfile $setting/external/external.X --gwas_thresh 0.01 --prune_thresh 0.1 --temp_basename ${exposure} --out ${outdir}${exposure}.prune

python3 scripts/make_mrtwin_infiles.py --trio_file $setting/trios/trio.txt --path_to_pheno $setting/trios/ --ext_bfile $setting/external/external.X --child_bfile $setting/trios/trio.child --p1_bfile $setting/trios/trio.father --p2_bfile $setting/trios/trio.mother --exposure_assoc ${outdir}${exposure}.gwas.PHENO1.glm.linear --outcome_assoc ${outdir}${outcome}.gwas.PHENO1.glm.linear --exposure_name ${exposure} --outcome_name ${outcome} --outname ${outdir}${pairname} --snps ${outdir}${exposure}.prune.in

python scripts/run_mrtwin_real_data.py --basename ${outdir}${pairname}
#
