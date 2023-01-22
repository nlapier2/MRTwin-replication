
## This script follows exactly the same simulation procedure as the function.R. 
# The only difference is to save the file instead to data path for downstream analysis.

source ("functions.R")
opts = commandArgs(trailingOnly=TRUE)
# opts=unlist(strsplit(commandArgs (TRUE)," "))
print(opts)
# library('rsvd')


if (length(opts) >= 1){
    seed = as.integer(opts[1])
} else {
    seed = 1
}
# Number of simulations
if (length(opts) >= 2){
    nsims = as.integer(opts[2])
    m = as.integer(opts[3])
    mcausal = as.integer(opts[4])
    n = as.integer(opts[5])
    # h2 = 0.1
    h2 = as.numeric(opts[6])
    se = as.numeric(opts[7])
    so = as.numeric(opts[8])
    gamma.ue = as.numeric(opts[9])
    gamma.uo = as.numeric(opts[10])
    ce = as.numeric(opts[11])
    ext.n = as.integer(opts[12])
    # ext.n = 100000
    ext_model = as.integer(opts[13])
    trio_model = as.integer(opts[14])
    ext.fst = as.numeric(opts[15])
    fst = as.numeric(opts[16])
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


ID=as.integer(opts[18])
cat ("*******", opts[18])
Threshold_flag = as.integer(opts[19]) # Flag: whether to use p-value threshold for feature selection 
if (Threshold_flag==1) {
   Threshold=0.05/m
} else {
   Threshold=1
}

cat("Threshold_flag is ", Threshold_flag, " threshold is: ", Threshold)
oFile=opts[20]

snpfilter="GWAS"
Nconf=0
Ncov=0
numtwins = 100
gwas_model = 0

pval = rep(-1,nsims)
pval2 = rep(-1,nsims)
pval3 = rep(-1,nsims)
pval4 = rep(-1,nsims)
pval5 = rep(-1,nsims)
pval6 = rep(-1,nsims)
pval7 = rep(-1,nsims)
pval8 = rep(-1,nsims)
pval9 = rep(-1,nsims)

# Relation between allele frequencies in trio and external
# 1: trio allele frequencies are equally related
# 2: trio allele frequencies are more closely related to one population than another
trio_ext =  1
# h2=0.01
set.seed (seed)
# browser()

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
cat ("# snpfilter:", snpfilter,"\n")



for (s in 1:nsims){
	true_betas = generate_true_beta (h2, m, mcausal)
    # true_betas = generate_true_beta_strong(m,coef=10)
	causal_status = true_betas$causal_status
	b = true_betas$b # strong effect
    # print(b)

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
    pca = generate_PC(ext.x)
    if (Ncov == 0){
        pcs=NULL
    }else{
         pcs = pca$x[,1:Ncov]
    }
   
    # pca = prcomp(ext.x,scale=TRUE,center = TRUE)
    # pcs = pca$x[,1:10]
    # pcs=rpca((ext.x),k=20)$x

    if (gwas_model == 0) {
        ext.betahat_e = run_gwas (ext.x, ext.e, covariates=pcs, alpha=Threshold, F=0,eflat=TRUE)
        ext.betahat_o = run_gwas (ext.x, ext.o, covariates=pcs, alpha=Threshold)
    } else if (gwas_model == 1) {
        ext.betahat_e = run_gwas (ext.x, ext.e, covariates = ext.poplabels, alpha=Threshold,F=0,eflat=TRUE)
        ext.betahat_o = run_gwas (ext.x, ext.o, covariates = ext.poplabels, alpha=Threshold)
    }

    if (trio_model == 0 ) {
        # Phased parental haplotypes (numtrios X 4)
        # D
        # if (trio_ext == 1) {
            # Trios are equally closely related to either population in reference
        x = rbinom (4*n * m, size = 1, prob = f)
        # } 
        # else if (trio_ext == 2) {
        #     # Trios are more closely related to one of the populations in reference
        #     x = rbinom (4*n * m, size = 1, prob = f1)
        # }

        x = matrix (x, 4*n, m, byrow = TRUE)

        # Population stratification will affect A,B,C and possibly D if PCs are used.
        # MR - estimate from external data
	if (ce==0){ # we wish to test the noisy version of external dataset to the calibration only
            ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter, causal_status)
            #ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o)
	}
	else{
	    ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter, causal_status)
        }
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
	if (ce==0){
            ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter)
	    #ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o)
	}
	else{
	    ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter)
	}


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
    x_orig = generate_offspring (xm1, xm2, xf1, xf2, numreps = 2)
    x = (x_orig[,,1])
    twinx=(x_orig[,,2])

    # Expsoure-outcome on offspring (we don't need exposure for our test)
    eo = generate_eo (x, b, se, so, causal_effect = ce, causal_status)
    e = eo$e
    o = eo$o 

    eo2 = generate_eo (twinx, b, se, so, causal_effect = ce, causal_status)
    e2 = eo2$e
    o2 = eo2$o 

    hits=ext.betahat_e$hits
    filters=ext.betahat_e$filters

    if (sum(filters&hits) != m){
        # print(filters&hits)
        print(paste0('WI1 detected',sum(filters&hits)))
        # next
    }

     if (trio_model == 1) {
        e = e + poplabels * gamma.ue
        o = o + poplabels * gamma.uo
    }
    pcs_trio = generate_PC(list(data=x,rotation=pca$rotation),type="transform")

    if (Ncov == 0){
        pcs_cov_trio = NULL
    }
    else{
        pcs_cov_trio = pcs_trio[,1:Ncov]
    }
    # print(pcs_cov_trio)
    betahat_e = run_gwas (x, e, covariates=pcs_cov_trio, alpha=1, F=0,eflat=TRUE) # no pval filter on the trio level
    hits2=betahat_e$hits
    filters2=betahat_e$filters
    if (sum(filters2&hits2) != m){
        # print(filters2&hits2)
        print(paste0('WI2 detected',sum(filters2&hits2)))
        # next
    }

    
    # if (ce==0){
    if (snpfilter=="ALL"){
        x_brump = x[,filters&filters2]
	print(paste0("X_brump shape is: ", dim(x_brump)))
        print(paste0("X shape is: ", dim(x)))
	brump_m = dim(x_brump)[2]
        xf = (xf1+xf2)
        xm = (xm1+xm2)
        # if (ce==0){
        xf_brump=xf[,filters&filters2]
        xm_brump=xm[,filters&filters2]
        # xf=xf[,filters]
    }
    else if (snpfilter=="GWAS"){
        
        x_brump = x[,filters&hits&filters2&hits2]
        # print(sum(filters&hits&filters2&hits2))
        if (sum(filters&hits&filters2&hits2)==0){
            next
        } else if (sum(filters&hits&filters2&hits2)==1){
            brump_m=1
        }
        else{
            brump_m = dim(x_brump)[2]
        }
        
        # print(paste("x_brump shape is", dim(x_brump)))
        # print(paste("x_brump shape is", (brump_m)))
        xf = (xf1+xf2)
        xm = (xm1+xm2)
        # if (ce==0){
        xf_brump=xf[,filters&hits&filters2&hits2]
        xm_brump=xm[,filters&hits&filters2&hits2]
        x = x[,hits&hits2]
        xf = xf[,hits&hits2]
        xm = xm[,hits&hits2]
        mrtwin_m = dim(x)[2]
        xm1 = xm1[,hits&hits2]
        xf1 = xf1[,hits&hits2]
        xm2 = xm2[,hits&hits2]
        xf2 = xf2[,hits&hits2]
    }
    else{
        print(paste0("snpfilter: ",snpfilter, " not implemented"))
        quit()
    }

    
   
   
    
    
    
    
    # }
	

    # Brunpton method
    
    # }
    pheno <- list(sibs1 = data.frame(x = e, y = o))
    fam_filter <- list(sibs1 = data.frame(x_brump), dads = data.frame(xf_brump), mums = data.frame(xm_brump))
    fam <- list(sibs1 = data.frame(x), dads = data.frame(xf), mums = data.frame(xm))
    # if (ce==0){# regress only legth(hits) number of snps)
    # print(paste0("Brumpton m is ",brump_m))
    #     b_sum_filter = brumpton_mr_trio(fam_filter, pheno,brump_m)
    #     b_sum = brumpton_mr_trio(fam, pheno,mrtwin_m)
    # }
    # else{
    #     b_sum = brumpton_mr_trio(fam, pheno,m)
    # }
    # Test statistic
    # if (ce==0){
    if (snpfilter=="GWAS"){
        ext.betahat_e$beta = ext.betahat_e$beta[hits&hits2]
        ext.betahat_e$se = ext.betahat_e$se[hits&hits2]
        ext.betahat_o$beta = ext.betahat_o$beta[hits&hits2]
        ext.betahat_o$se = ext.betahat_o$se[hits&hits2]
    }
    else{
        ext.betahat_e$beta = ext.betahat_e$beta[hits]
        ext.betahat_e$se = ext.betahat_e$se[hits]
        ext.betahat_o$beta = ext.betahat_o$beta[hits]
        ext.betahat_o$se = ext.betahat_o$se[hits]
    }

    betahatEO=rbind(ext.betahat_e$beta,ext.betahat_e$se,ext.betahat_o$beta, ext.betahat_o$se)
    eo1 = rbind(c(e),c(o))
    eo2 = rbind(c(e2),c(o2))
    cat ("writing betaEO file", dim(betahatEO), "\n")
    write.table(betahatEO,paste0(oFile,"/",ID,".",seed, ".", s,".betahatEO"),row.names=FALSE, col.names=FALSE)
    cat("writing p1 file",dim(xf),"\n")
    write.table(xf,paste0(oFile,"/",ID,".",seed, ".", s,".pa1geno"),row.names=FALSE, col.names=FALSE)
    cat("writing p2 file",dim(xm),"\n")
    write.table(xm,paste0(oFile,"/",ID,".", seed, ".", s,".pa2geno"),row.names=FALSE, col.names=FALSE)
    cat("writing x file",dim(x),"\n")
    write.table(x,paste0(oFile,"/",ID,".",seed, ".", s,".childgeno"),row.names=FALSE, col.names=FALSE)
    cat("writing twinx file",dim(twinx[,hits&hits2]),"\n")
    write.table(twinx[,hits&hits2],paste0(oFile,"/",ID,".",seed, ".", s,".child2geno"),row.names=FALSE, col.names=FALSE)
    cat("writing exposure outcome trio1 file",dim(eo1),"\n")
    write.table(eo1,paste0(oFile,"/",ID,".",seed, ".", s,".childpheno"),row.names=FALSE, col.names=FALSE)
    cat("writing exposure outcome trio2 file",dim(eo2),"\n")
    write.table(eo2,paste0(oFile,"/",ID,".", seed, ".", s,".child2pheno"),row.names=FALSE, col.names=FALSE)
    
}



