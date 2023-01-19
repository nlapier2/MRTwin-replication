library(MendelianRandomization)

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
corfile = args[2]

cormat = as.matrix(read.table(corfile))
betahatEO = as.matrix(read.table(infile))
beta_e = betahatEO[1,]
se_e = betahatEO[2,]
beta_o = betahatEO[3,]
se_o = betahatEO[4,]

mr_input_obj = MendelianRandomization::mr_input(bx = beta_e, bxse= se_e, by = beta_o, byse = se_o, correlation = cormat)
ivw_results = MendelianRandomization::mr_ivw(mr_input_obj, correl = TRUE)
egger_results = MendelianRandomization::mr_egger(mr_input_obj, correl = TRUE)
median_results = MendelianRandomization::mr_median(mr_input_obj)
ivw_pval = ivw_results$Pvalue
egger_pval = egger_results$Causal.pval
median_pval = median_results$Pvalue

cat(ivw_pval, '\n')
cat(egger_pval, '\n')
cat(median_pval, '\n')
#

