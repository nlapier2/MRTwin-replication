#!/usr/bin/env Rscript

library("stringr")
# helper function for brumpton method
mr_ivw <- function(beta_e, se_e, beta_o, se_o)
{
    num = (beta_e * beta_o) %*% (1/(se_o^2))
    denom = (beta_e^2) %*% (1/(se_o^2))
    mr_est = num / denom
    mr_se = sqrt(1 / denom)
    mr_stat = mr_est^2 / (mr_se^2)
    mr_pval = 1 - pchisq(mr_stat, 1)
    return(list("eo" = mr_est, "se" = mr_se, "t" = mr_stat, "p" = mr_pval))
}


# brumpton trio method
brumpton_mr_trio <- function(fam, phen, nsnp, names, threshold)
{
    bx <- rep(NA,nsnp)
    by <- rep(NA,nsnp)
    sex <- rep(NA,nsnp)
    sey <- rep(NA,nsnp)
    for(i in 1:nsnp)
    {
        #sdiffgx <- fam$sibs1[,i] * bx[i] - fam$sibs2[,i] * bx[i]
        #sdiffx <- phen$sibs1$x - phen$sibs2$x
        if(names[1] %in% c("derived_chd", "derived_Stroke", "diabetes", "asthma"))
{
           mod <- summary(glm(as.factor(phen$sibs1$x) ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i], family = "binomial"))$coefficients
 
        }
        else
{
	modsum <- summary(lm(phen$sibs1$x ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))
	#print(modsum)
	mod <- modsum$coefficients
	modpval = mod[2,4]
	f = as.numeric(modsum$f['value'])
        #mod <- summary(lm(phen$sibs1$x ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))$coefficients
} 
   if (nrow(mod) ==  1 || modpval > threshold) {  # f < threshold) { 
        bx[i] <- 0
        sex[i] <- 1e6
    }
    else {
        bx[i] <- mod[2,1]
            sex[i] <- mod[2,2]
        }
        if(names[2] %in% c("derived_chd", "derived_Stroke", "diabetes", "asthma"))
        {
           mod <- summary(glm(as.factor(phen$sibs1$y) ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i], family = "binomial"))$coefficients

        }
else
{

  mod <- summary(lm(phen$sibs1$y ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))$coefficients

}
    #print(mod)
    if (nrow(mod) ==  1 || modpval > threshold) {  #f < threshold) {
        by[i] <- 0
        sey[i] <- 1e6
    }
    else {
        by[i] <- mod[2,1]
            sey[i] <- mod[2,2]
    }
    }
    if (min(bx) == 0 && max(bx) == 0) {
        return(list('p'=1.0))
    }
    return(mr_ivw(bx,sex, by, sey))
}

# read in command line args
args = commandArgs(trailingOnly=TRUE)
basename = args[1]
names = unlist(str_split(args[3], "~"))
# run brumpton
if (!(file.exists(paste0(basename, '.childgeno')))) {
	quit()
}
betahatEO <- read.table(paste0(basename, '.betahatEO'))
num_snps = length(betahatEO) - 1
threshold = min(1 / num_snps, 0.05)  # 10
childgeno <- read.table(paste0(basename, '.childgeno'))
pa1geno   <- read.table(paste0(basename, '.pa1geno'))
pa2geno   <- read.table(paste0(basename, '.pa2geno'))
childEO   <- read.table(paste0(basename, '.childpheno'))
e <- head(childEO, 1)
o <- tail(childEO, 1)
pheno <- list(sibs1 = data.frame(x = as.numeric(e), y = as.numeric(o)))
fam <- list(sibs1=childgeno, dads=pa1geno, mums=pa2geno)
pval = brumpton_mr_trio(fam, pheno, num_snps, names, threshold)$p
#print(warnings())
cat(pval,'\n')
#
