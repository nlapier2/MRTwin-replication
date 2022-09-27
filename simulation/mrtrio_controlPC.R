source ("functions.R")
# opts = unlist(strsplit("1.00 500.00 200.00 100.00 1000.00 0.0001 0.9999 1 0.1 0.1 0.0 5000.00 1.00 1.00 0.10 0.10", " "))
opts=unlist(strsplit(commandArgs (TRUE)," "))
print(opts)
if (length(opts) >= 1){
    seed = as.integer(opts[1])
} else {
    seed = 1
}
# Number of simulations
if (length(opts) >= 2){
    nsims = as.integer(opts[2]) # number of total simulation
    m = as.integer(opts[3]) # number of total snps
    mcausal = as.integer(opts[4]) # number of effective snps 
    n = as.integer(opts[5]) # trio sample size
    h2 = as.numeric(opts[6]) # heritabiilty from genetics to exposure (Instrumental strength)
    se = as.numeric(opts[7]) # noise level on exposure trait
    so = as.numeric(opts[8]) # noise level on outcome trait
    gamma.ue = as.numeric(opts[9])
    gamma.uo = as.numeric(opts[10])
    ce = as.numeric(opts[11]) # causal effect from exposure to outcome
    ext.n = as.integer(opts[12]) # external sample size
    ext_model = as.integer(opts[13]) # external model = 1 will generate two population groups with equal size; external model = 0 will generate single population
    trio_model = as.integer(opts[14]) # similar to ext_model, suggest to set with the same value as ext_model
    ext.fst = as.numeric(opts[15]) # Fst in external population (ext_model=1)
    fst = as.numeric(opts[16]) # Fst in trio population (trio_model=1)
} else {
    nsims = 1000
    # Number of SNPs
    m = 1

    # Number of causal SNPs (smallest is 2. if we try 1, R converts matrices to vectors)
    mcausal = 1

    # Number of trios
    n = 100

    # Number of digital twins

    # Number of external genotypes
    ext.n = 1000

    # Genetic variance of all SNPs on e
    h2 = 0.5
    # Environmental variance of e
    se = 0.5
    # Environmental variance of o
    so = 0.5


    # Effect of population stratification on exposure 
    # Relevant only under structured models
    gamma.ue = 1
    #
    # Effect of population stratification on outcome 
    # Relevant only under structured models
    gamma.uo = 1

    # Causal effect
    ce = 0.5

    # Model for the external reference data (simple : no structure/ structured: has population structure)
    #ext_model = "simple"
    # 1 means structured, 0 means simple
    ext_model = 0

    # Model for the effect size estimates/GWAS estimates (simple: no correction for population structure/ oracle: correction using true labels / pc: correction using PCs )
    
    #gwas_model = "oracle" => 1
    #gwas_model = "simple" => 0

    #trio_model = "simple"
    # 1 means structured, 0 means simple
    trio_model = 0

    ext.fst = 0.05
    fst = 0.05
}


numtwins = 100
gwas_model = 0

pval = rep(NA,nsims)
pval2 = rep(NA,nsims)
pval3 = rep(NA,nsims)

# Relation between allele frequencies in trio and external
# 1: trio allele frequencies are equally related
# 2: trio allele frequencies are more closely related to one population than another
trio_ext =  1

set.seed (seed)

cat ("# seed:", seed,"\n")
cat ("# h2:", h2,"\n")
cat ("# m:", m,"\n")
cat ("# mcausal:", mcausal,"\n")
cat ("# num_trios:", n,"\n")
cat ("# num_twins:", numtwins,"\n")
cat ("# num_ext:", ext.n, "\n")
cat ("# ext_model:", ext_model,"\n")
cat ("# gwas_model:", gwas_model,"\n")
cat ("# trio_model:", trio_model,"\n")
cat ("# trio_ext:", trio_ext,"\n")
cat ("# causal_effect:", ce,"\n")
cat ("# effect_of_population_on_exposure:", gamma.ue,"\n")
cat ("# effect_of_population_on_outcome:", gamma.uo,"\n")
cat ("# exposure noise:", se,"\n")
cat ("# outcome noise:", so,"\n")

snpfilter="gwas"

for (s in 1:nsims){
	true_betas = generate_true_beta (h2, m, mcausal)
	causal_status = true_betas$causal_status
	b = true_betas$b
	f = generate_af (m, "uniform")


    if (ext_model == 0){
        # Reference genotypes
        # A
        ext.x = rbinom (ext.n*m, size = 2, prob = f) 
        ext.x = (matrix (ext.x, ext.n, m, byrow = TRUE))

        # B
        eo_ext = generate_eo (ext.x, b, se, so, causal_effect = ce, causal_status= causal_status)
        ext.e = eo_ext$e
        ext.o = eo_ext$o
    } else if (ext_model == 1) {

        

        cat ("# fst:", ext.fst, "\n")
        ext.f = generate_af (m, "uniform")
        ext.f1 = rbeta (m, ext.f*(1-ext.fst)/ext.fst,(1-ext.f)*(1-ext.fst)/ext.fst) 
        ext.f2 = rbeta (m, ext.f*(1-ext.fst)/ext.fst,(1-ext.f)*(1-ext.fst)/ext.fst) 

        rfst = get_fst (ext.f1,ext.f2)
        cat ("# Realized fst:", rfst, "\n")
        #    f1 = generate_af (m, "uniform")
        #    f2 = generate_af (m, "uniform")
        ext.poplabels = rep(1, ext.n) 
        ext.poplabels[1:round(0.5*ext.n)] = 0
        ext.n1 = length(which(ext.poplabels==0))
        ext.n2 = length(which(ext.poplabels==1))

        # Reference genotypes
        # A
        ext.x1 = rbinom (ext.n1*m, size = 2, prob = ext.f1) 
        ext.x2 = rbinom (ext.n2*m, size = 2, prob = ext.f2) 
        ext.x = c(ext.x1,ext.x2)
        ext.x = (matrix (ext.x, ext.n, m, byrow = TRUE))

        # B
        eo_ext = generate_eo (ext.x, b, se, so, causal_effect = ce, causal_status=causal_status)
        ext.e = eo_ext$e + ext.poplabels * gamma.ue
        ext.o = eo_ext$o + ext.poplabels * gamma.uo

    }

    # July 9th 2022: remove scale for matrix (this could cause issue for digital twin)
    # C
    pca = prcomp(ext.x,scale=TRUE,center = TRUE)
    pcs = pca$x[,1:5]

    if (gwas_model == 0) {
        ext.betahat_e = run_gwas (ext.x, ext.e, covariates=pcs)

        ext.betahat_o = run_gwas (ext.x, ext.o, covariates=pcs)
    } else if (gwas_model == 1) {
        ext.betahat_e = run_gwas (ext.x, ext.e, covariates = ext.poplabels,covariates=pcs)
        ext.betahat_o = run_gwas (ext.x, ext.o, covariates = ext.poplabels,covariates=pcs)
    }

    if (trio_model == 0 ) {
        # Phased parental haplotypes (numtrios X 4)
        # D
        if (trio_ext == 1) {
            # Trios are equally closely related to either population in reference
            x = rbinom (4*n * m, size = 1, prob = f)
        } else if (trio_ext == 2) {
            # Trios are more closely related to one of the populations in reference
            x = rbinom (4*n * m, size = 1, prob = f1)
        }

        x = matrix (x, 4*n, m, byrow = TRUE)

        # Population stratification will affect A,B,C and possibly D if PCs are used.
        # MR - estimate from external data

	    ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter, causal_status)

	# Maternal haplotypes
        xm1 = x[1:n,]
        xm2 = x[(n+1):(2*n),]

        # Paternal haplotypes
        xf1 = x[(2*n+1):(3*n),]
        xf2 = x[(3*n+1):(4*n),]
    } else if (trio_model == 1) {
        

        cat ("# fst:", fst, "\n")
        f = generate_af (m, "uniform")
        # f1 = rbeta (m, f*(1-fst)/fst,(1-f)*(1-fst)/fst) 
        # f2 = rbeta (m, f*(1-fst)/fst,(1-f)*(1-fst)/fst) 
        # July 9th 2022: make the maf of external and trio consistent
        f1 = ext.f1
        f2 = ext.f2

        rfst = get_fst (f1,f2)
        cat ("# Realized fst:", rfst, "\n")

        # Phased parental haplotypes (numtrios X 4)
        # D
        x1 = rbinom (2*n * m, size = 1, prob = f1)
        x1 = matrix (x1, 2*n, m, byrow = TRUE)
        x2 = rbinom (2*n * m, size = 1, prob = f2)
        x2 = matrix (x2, 2*n, m, byrow = TRUE)
        # Population stratification will affect A,B,C and possibly D if PCs are used.
        # MR - estimate from external data

	    ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter)
        n0 = round(0.5*n)

        poplabels = rep(1, n) 
        poplabels[1:round(n0)] = 0
        n1 = length(which(poplabels==0))
        n2 = length(which(poplabels==1))
            
        # Maternal haplotypes
        xm1 = rbind(x1[1:n0,],x2[1:n0,])
        xm2 = rbind(x1[(n0+1):n,],x2[(n0+1):n,])

        # Paternal haplotypes
        xf1 = rbind(x1[(n+1):(n+n0),],x2[(n+1):(n+n0),])
        xf2 = rbind(x1[(n+n0+1):(2*n),],x2[(n+n0+1):(2*n),])
    }

    #for (s in 1:nsims){
    # Initial individuals
    x = generate_offspring (xm1, xm2, xf1, xf2, numreps = 1)
    x = (x[,,1])

    # Expsoure-outcome on offspring (we don't need exposure for our test)
    eo = generate_eo (x, b, se, so, causal_effect = ce, causal_status)
    e = eo$e
    o = eo$o 
    hits=ext.betahat_e$hits
    if (ce==0){
        x = x[,hits]
    }
	if (trio_model == 1) {
		e = e + poplabels * gamma.ue
		o = o + poplabels * gamma.uo
	}

    # Brunpton method
    xf = (xf1+xf2)
    xm = (xm1+xm2)
    if (ce==0){
        xf=xf[,hits]
	    xm=xm[,hits]
    }
    pheno <- list(sibs1 = data.frame(x = e, y = o))
    fam <- list(sibs1 = data.frame(x), dads = data.frame(xf), mums = data.frame(xm))

    if (ce==0){# regress only legth(hits) number of snps)
        b_sum = brumpton_mr_trio(fam, pheno,length(hits))
    }
    else{
        b_sum = brumpton_mr_trio(fam, pheno,m)
    }
    t2 = b_sum$t
    p2 = b_sum$p

    # Test statistic
    if (ce==0){
        t = get_stat (x, ext.betahat_e, ext.betahat_eo, o, filter=TRUE)
    }
    else{
        t = get_stat (x, ext.betahat_e, ext.betahat_eo, o)
    }
    # Digital twins
    # twins = generate_offspring (xm1, xm2, xf1, xf2, numreps = numtwins)
    twint = rep(NA, numtwins)
    twint2 = rep(NA, numtwins)
    for ( i in 1:numtwins) {
        twins = generate_offspring (xm1, xm2, xf1, xf2, numreps = 1)
        twinx = (twins[,,1])
	if(ce==0){
	    twinx = twinx[,hits]
	    twint[i] = get_stat (twinx, ext.betahat_e, ext.betahat_eo, o, filter=TRUE)
	}
	else{
            twint[i] = get_stat (twinx, ext.betahat_e, ext.betahat_eo, o)
	}
    }
    
    mr_pval = mr_ivw(ext.betahat_e$beta,ext.betahat_e$se,ext.betahat_o$beta,ext.betahat_o$se)$p
    pval[s] = (1 + length(which (twint > t)))/(1+numtwins)
    pval2[s] = p2
#    pval2[s] = (1 + length(which (twint2 > rep(t2,length(twint2)))))/(1+numtwins)
    pval3[s] = mr_pval

    cat ("Replicate: ", s, "mrtrio pval = ", pval[s],"\n") 
    cat ("Replicate: ", s, "brumpton pval  = ", pval2[s],"\n") 
    cat ("Replicate: ", s, "ivw pval  = ", pval3[s],"\n") 
}

#Test calibration
cal_05=get_calibration(pval, 0.05)
cal2_05=get_calibration(pval2, 0.05)
cal3_05=get_calibration(pval3, 0.05)
cat ("## mrtrio Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal_05$n0, cal_05$n, cal_05$p,cal_05$pp,"\n")

cat ("## Brumpton Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal2_05$n0, cal2_05$n, cal2_05$p, cal2_05$pp,"\n")

cat ("## IVW Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal3_05$n0, cal3_05$n, cal3_05$p, cal3_05$pp,"\n")

