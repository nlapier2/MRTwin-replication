
library('MASS')
library(MendelianRandomization)
# get summary statistics when there are only two populations involved
gwas_summary = function(N1, N2, f2, se_e, se_o, ce, beta_e, beta_o, gamma.ue, gamma.uo){
    
    mean_be = beta_e + (2*N2/(N1+N2))*gamma.ue*f2 
    mean_bo = beta_o + ce*beta_e + (2*N2/(N1+N2))*(ce*gamma.ue+gamma.uo)*f2
    m = length(mean_be)
    c11 = (1/(N1+N2))*se_e*diag(m)
    c12 = ce*(1/(N1+N2))*se_e*diag(m)
    c22 = (1/(N1+N2))*(ce**2*se_e+se_o)*diag(m)
    beta_hat = mvrnorm(n=1,c(mean_be, mean_bo),  rbind(cbind(c11,c12), cbind(c12,c22)))
    beta_hat_e = beta_hat[1:m]
    beta_hat_o = beta_hat[(m+1):(2*m)]
    # browser()
    betahat_e = list(beta=beta_hat_e, se = sqrt(diag(c11)))
    betahat_o = list(beta=beta_hat_o, se = sqrt(diag(c22)))
    return (list (betahat_e=betahat_e, betahat_o=betahat_o))
}



get_fst = function (p1, p2) {	
	n = sum( (p1-p2)^2)
	d = sum( p1*(1-p2)+p2*(1-p1))
	f = 0
	if (d>0){
		f = n/d
	}
	return (f)
}




get_calibration = function (pval, alpha){
    n = length(pval)
    p = length(which(pval<alpha))/n
    se = sqrt(alpha*(1-alpha)/n)
    ose = sqrt(p*(1-p)/n)
    z = (p-alpha)/se
    pp = 2 * pnorm(-abs(z))

    return (list (n = n, n0 = length( which (pval < alpha) ) , p=p, se=se, ose=ose, z=z,pp=pp) )
}

generate_offspring = function (xm1, xm2, xf1, xf2, numreps = 100){
    n = nrow (xm1)
    m = ncol (xm1)
    off = array (NA, dim = c(n, m, numreps))
    for ( i in 1:numreps){
        tm = rbinom(n*m, size = 1, prob=0.5)
        tm = matrix(tm, n, m)
        tf = rbinom(n*m, size = 1, prob=0.5)
        tf = matrix(tf, n, m)

        xf = xf1*(1-tf) + xf2 * tf
        xm = xm1*(1-tm) + xm2 * tm

        off[,,i] = xf + xm
    }
    return (off)
}

generate_eo = function (x, b, se, so, causal_effect, causal_status) {
    e = x[,which(causal_status==1)]%*%b[which(causal_status==1)] + rnorm (nrow(x), mean = 0, sd = sqrt(se))
    o = e*causal_effect + rnorm (nrow(x), mean = 0, sd = sqrt(so))
    return (list(e=e, o=o))
}

run_gwas = function (x, pheno, covariates = NULL, debug = FALSE,alpha=0.05) {
    beta = rep(NA, ncol(x))
    se = rep(NA, ncol(x))
    p = rep(NA, ncol(x))
    M = dim(x)[2]
    # test increasing the threshold
    thresh = alpha
#    thresh = alpha/M
    hits = c()
    for (i in 1:(ncol(x))){
        if (is.null(covariates)){ 
            l = lm (pheno~x[,i])
        } else {
            l = lm (pheno~x[,i] + covariates )
        }
        if (debug == TRUE) {
            print(summary(l))
        }
        if (nrow(summary(l)$coefficients) == 1) { 
            # Ideally, we will filter this SNP out
            # Instead setting p = 1 and beta = 0
            p[i] = 1
            beta[i] = 0
            se[i] = 1e6
        } else { 
            beta[i] = summary(l)$coefficients[2,1]
            se[i] = summary(l)$coefficients[2,2]
            p[i] = summary(l)$coefficients[2,4]
        }
	if (p[i] <= thresh){
	    hits <- c(hits,i)
	}
    }
#    p1 = which (p<0)
#    beta = beta[p1]
#    se = se [p1]
#    p = p[p1]
    return (list (beta=beta, se=se, p=p, hits=hits))
}

generate_af = function (m, distribution = "same"){
    if (distribution == "same") {
        f = rep (0.25, m)
    } else {    
        f = runif (m, min = 0.1, max = 0.5)
    } 
    return (f)    
}

generate_true_beta = function (h2, m, mcausal) { 
    causal_status = rep(0,m)
    causal_status[1:mcausal] = 1

    b = rnorm (m, mean = 0, sd = sqrt(h2/mcausal))
    b = b * causal_status

    return (list ( causal_status = causal_status, b = b))
}


# snpfilter : which set of snps to use (oracle: use true causals, all: use all)
get_mr = function (betahat_e, betahat_o, snpfilter = "gwas", causal_status = NULL, weight='ivw') {
   # Ideally, we would use IVW or the median estimator but it should not matter too much for our purpose
   if (weight == 'uniform'){
       if (snpfilter == "all"){  
       betahat_eo = mean(betahat_o$beta/betahat_e$beta)
        } else if (snpfilter == "oracle") {
            whichsnp = betahat_e$p < 1 & betahat_o$p < 1
            betahat_eo = betahat_o$beta/betahat_e$beta
            betahat_eo = betahat_eo[which(causal_status==1) & whichsnp]
            betahat_eo = mean(betahat_eo)
        }
   }
   else if (weight == 'ivw'){
       if (snpfilter=='gwas'){
          hits = betahat_e$hits
	  betahat_e$beta=betahat_e$beta[hits]
	  betahat_e$se=betahat_e$se[hits]
	  betahat_o$beta=betahat_o$beta[hits]
	  betahat_o$se=betahat_o$se[hits]
       }
       betahat_eo = mr_ivw(betahat_e$beta, betahat_e$se, betahat_o$beta, betahat_o$se)$eo[1]
   }
   return (betahat_eo)
}




get_stat = function (x, betahat_e,  betahat_eo, o, filter=FALSE){
    beta = betahat_e$beta
    if (filter){
        hits = betahat_e$hits
        beta = beta[hits]
    }
    if(is.null(nrow(x)) || ncol(x)==1){ # 1 snp case
       t = - sum((o - (x * betahat_e$beta ) * betahat_eo)^2)
       return (t)
    }
    t = - sum((o - (x %*% beta ) * betahat_eo)^2)
    return (t)
}

get_pc_scores = function (x, k = 1) { 
    res.pca = prcomp (x , scale = TRUE)
    res.ind = get_pca_ind (res.pca)
    return (res.ind$coord[,1:k])
}

qqunif = function (p) { 
    qqplot (p, runif(length(p)))
}


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

mr_egger<- function(beta_e, se_e, beta_o, se_o)
{
  mr_input_obj = MendelianRandomization::mr_input(bx = beta_e, bxse= se_e, by = beta_o, byse = se_o)
  egger_results = MendelianRandomization::mr_egger(mr_input_obj)
  mr_est = egger_results$Estimate
  mr_se = egger_results$StdError.Est
  p_val = egger_results$Causal.pval
  return(list("eo" = mr_est, "se" = mr_se, "p" = p_val))
}

#fam needs to have sibs1, dads and mums
#pheno x and y from sibs
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
	options(error = recover)
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

