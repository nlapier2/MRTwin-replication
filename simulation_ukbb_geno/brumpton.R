#!/usr/bin/env Rscript


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
brumpton_mr_trio <- function(fam, phen, nsnp)
{
    bx <- rep(NA,nsnp)
    by <- rep(NA,nsnp)
    sex <- rep(NA,nsnp)
    sey <- rep(NA,nsnp)
    for(i in 1:nsnp)
    {
        #sdiffgx <- fam$sibs1[,i] * bx[i] - fam$sibs2[,i] * bx[i]
        #sdiffx <- phen$sibs1$x - phen$sibs2$x
        mod <- summary(lm(phen$sibs1$x ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))$coefficients
    if (nrow(mod) ==  1){
        bx[i] <- 0
        sex[i] <- 1e6
    }
    else {
        bx[i] <- mod[2,1]
            sex[i] <- mod[2,2]
        }
    mod <- summary(lm(phen$sibs1$y ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))$coefficients
    if (nrow(mod) ==  1){
        by[i] <- 0
        sey[i] <- 1e6
    }
    else {
        by[i] <- mod[2,1]
            sey[i] <- mod[2,2]
    }
    }
    return(mr_ivw(bx,sex, by, sey))
}


# read in command line args
args = commandArgs(trailingOnly=TRUE)
basename = args[1]
num_sims = as.numeric(args[2])
threshold = 0.05
positives = 0
successful_sims = 0
# run brumpton
for (i in 1:num_sims) {
	sim=i-1
	if (!(file.exists(paste0(basename, sim, '.child1geno')))) {
		next
	}
	successful_sims = successful_sims + 1
        betahatEO <- read.table(paste0(basename, sim, '.betahatEO'))
        num_snps = length(betahatEO) - 1
	childgeno <- read.table(paste0(basename, sim, '.child1geno'))
	pa1geno   <- read.table(paste0(basename, sim, '.pa1geno'))
	pa2geno   <- read.table(paste0(basename, sim, '.pa2geno'))
	childEO   <- read.table(paste0(basename, sim, '.child1pheno'))
	e <- head(childEO, 1)
	o <- tail(childEO, 1)
	pheno <- list(sibs1 = data.frame(x = as.numeric(e), y = as.numeric(o)))
	fam <- list(sibs1=childgeno, dads=pa1geno, mums=pa2geno)
	pval = brumpton_mr_trio(fam, pheno, num_snps)$p
	cat(pval,'\n')
	if(pval < threshold) {
		positives = positives + 1
	}
}
positive_rate = positives / successful_sims # num_sims
cat(paste0('FPR/Power: ', positive_rate), '\n')
#
