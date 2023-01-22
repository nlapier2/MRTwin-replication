source ("functions.R")
opts = commandArgs(trailingOnly=TRUE)
print(opts)


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
    h2 = as.numeric(opts[6])
    se = as.numeric(opts[7])
    so = as.numeric(opts[8])
    gamma.ue = as.numeric(opts[9])
    gamma.uo = as.numeric(opts[10])
    ce = as.numeric(opts[11])
    ext.n = as.integer(opts[12])
    ext_model = as.integer(opts[13])
    trio_model = as.integer(opts[14])
    ext.fst = as.numeric(opts[15])
    fst = as.numeric(opts[16])
    snpfilter=opts[17]
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

pval10 = rep(-1,nsims)
pval11 = rep(-1,nsims)
pval12 = rep(-1,nsims)
pval13 = rep(-1,nsims)



weight="ivw"

# Relation between allele frequencies in trio and external
# 1: trio allele frequencies are equally related
# 2: trio allele frequencies are more closely related to one population than another
trio_ext =  1

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

# snpfilter="gwas"

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
        ext.betahat_e = run_gwas (ext.x, ext.e, covariates=pcs, alpha=0.05/m, F=0,eflat=TRUE)
        ext.betahat_o = run_gwas (ext.x, ext.o, covariates=pcs, alpha=0.05/m)
    } else if (gwas_model == 1) {
        ext.betahat_e = run_gwas (ext.x, ext.e, covariates = ext.poplabels, alpha=0.05/m,F=0,eflat=TRUE)
        ext.betahat_o = run_gwas (ext.x, ext.o, covariates = ext.poplabels, alpha=0.05/m)
    }
    mae_betae=mean((ext.betahat_e$beta-b))
    mae_betao=mean((ext.betahat_o$beta-0))
    selected=(ext.betahat_e$hits)&(ext.betahat_e$filters)
    smae_betae=mean((ext.betahat_e$beta[selected]-b[selected]))
    smae_betao=mean((ext.betahat_o$beta[selected]))
    corr_eo=cor(ext.betahat_e$beta, ext.betahat_o$beta, method = 'pearson')
    corr_eo_sel=cor(ext.betahat_e$beta[selected], ext.betahat_o$beta[selected], method = 'pearson')
    print(paste0("MAE beta_e est: ",mae_betae))
    print(paste0("MAE beta_o est: ",mae_betao))
    print(paste0("est beta_e beta_o corr: ",corr_eo))
    print(paste0("selected MAE beta_e est: ",smae_betae))
    print(paste0("selected MAE beta_o est: ",smae_betao))
    print(paste0("selected est beta_e beta_o corr: ",corr_eo_sel))


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
            ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter, causal_status=causal_status,weight=weight)
            #ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o)
	}
	else{
	    ext.betahat_eo = get_mr (ext.betahat_e, ext.betahat_o, snpfilter = snpfilter, causal_status=causal_status,weight=weight)
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
    x = generate_offspring (xm1, xm2, xf1, xf2, numreps = 1)
    x = (x[,,1])

    # Expsoure-outcome on offspring (we don't need exposure for our test)
    eo = generate_eo (x, b, se, so, causal_effect = ce, causal_status)
    e = eo$e
    o = eo$o 
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
    betahat_e = run_gwas (x, e, covariates=pcs_cov_trio, alpha=1, F=10,eflat=TRUE) # no pval filter on the trio level
    betahat_o = run_gwas (x, o, covariates=pcs_cov_trio, alpha=1, F=10,eflat=TRUE) # no pval filter on the trio level
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
    print(paste0("Brumpton m is ",brump_m))
    b_sum_filter = brumpton_mr_trio(fam_filter, pheno,brump_m)
    b_sum_tdist = brumpton_mr_trio_tdist(fam, pheno,mrtwin_m)
    b_sum = brumpton_mr_trio(fam, pheno,mrtwin_m)
    # }
    # else{
    #     b_sum = brumpton_mr_trio(fam, pheno,m)
    # }
    t2 = b_sum_filter$t
    p2 = b_sum_filter$p

    p7 = b_sum$p
    # p8 = b_sum_tdist$p
    # print("start calculating ivw pval")
    

    ext.betahat_e.beta.filter=ext.betahat_e$beta[hits&hits2&filters2&filters]
    ext.betahat_o.beta.filter=ext.betahat_o$beta[hits&hits2&filters2&filters]
    ext.betahat_e.se.filter=ext.betahat_e$se[hits&hits2&filters2&filters]
    ext.betahat_o.se.filter=ext.betahat_o$se[hits&hits2&filters2&filters]

    # Test statistic
    # if (ce==0){
    if (snpfilter=="GWAS"){
        ext.betahat_e$beta = ext.betahat_e$beta[hits]
        ext.betahat_e$se = ext.betahat_e$se[hits]
        ext.betahat_o$beta = ext.betahat_o$beta[hits]
        ext.betahat_o$se = ext.betahat_o$se[hits]
    }
    else{
        ext.betahat_e$beta = ext.betahat_e$beta[hits&filter]
        ext.betahat_e$se = ext.betahat_e$se[hits&filter]
        ext.betahat_o$beta = ext.betahat_o$beta[hits&filter]
        ext.betahat_o$se = ext.betahat_o$se[hits&filter]
    }

    ## replace IVW with uniform weight

    ## testing on trio set to avoid double dipping
    # mr_pval_t = mr_ivw_tdist(betahat_e$beta[hits],betahat_e$se[hits],betahat_o$beta[hits],betahat_o$se[hits])$p
    mr_pval = mr_ivw(betahat_e$beta[hits],betahat_e$se[hits],betahat_o$beta[hits],betahat_o$se[hits])$p

    ## newly added Nov 29th
    # print("start calculating egger pval")
    egger_pval = mr_egger(betahat_e$beta[hits],betahat_e$se[hits],betahat_o$beta[hits],betahat_o$se[hits])$p
    # print("start calculating median pval")
    median_pval = mr_median(betahat_e$beta[hits],betahat_e$se[hits],betahat_o$beta[hits],betahat_o$se[hits])$p
    # print("start calculating mode pval")
    mode_pval = mr_mode(betahat_e$beta[hits],betahat_e$se[hits],betahat_o$beta[hits],betahat_o$se[hits])$p


    mr_pval_filter = mr_ivw(ext.betahat_e.beta.filter,ext.betahat_e.se.filter,ext.betahat_o.beta.filter,ext.betahat_o.se.filter)$p
    ## newly added Nov 29th
    # print("start calculating egger pval")
    egger_pval_filter = mr_egger(ext.betahat_e.beta.filter,ext.betahat_e.se.filter,ext.betahat_o.beta.filter,ext.betahat_o.se.filter)$p
    # print("start calculating median pval")
    median_pval_filter = mr_median(ext.betahat_e.beta.filter,ext.betahat_e.se.filter,ext.betahat_o.beta.filter,ext.betahat_o.se.filter)$p
    # print("start calculating mode pval")
    mode_pval_filter = mr_mode(ext.betahat_e.beta.filter,ext.betahat_e.se.filter,ext.betahat_o.beta.filter,ext.betahat_o.se.filter)$p
    
    t = get_stat (x, ext.betahat_e, ext.betahat_eo, o)
    # }
    # else{
    #     t = get_stat (x, ext.betahat_e, ext.betahat_eo, o)
    # }
    # Digital twins
    # twins = generate_offspring (xm1, xm2, xf1, xf2, numreps = numtwins)
    twint = rep(NA, numtwins)
    twint2 = rep(NA, numtwins)
    for ( i in 1:numtwins) {
        twins = generate_offspring (xm1, xm2, xf1, xf2, numreps = 1)
        twinx = (twins[,,1])
	# if(ce==0){
	    # twinx = twinx[,hits]
	    twint[i] = get_stat (twinx, ext.betahat_e, ext.betahat_eo, o)
	# }
	# else{
    #     twint[i] = get_stat (twinx, ext.betahat_e, ext.betahat_eo, o,filter=FALSE)
	# }
    }
    
    
    pval[s] = (1 + length(which (twint > t)))/(1+numtwins)
    pval2[s] = p2
#    pval2[s] = (1 + length(which (twint2 > rep(t2,length(twint2)))))/(1+numtwins)
    pval3[s] = mr_pval
    pval4[s] = egger_pval
    pval5[s] = median_pval
    pval6[s] = mode_pval
    pval7[s] = p7
    pval10[s] = mr_pval_filter
    pval11[s] = egger_pval_filter
    pval12[s] = median_pval_filter
    pval13[s] = mode_pval_filter
    # pval8[s] = p8
    # pval9[s] = mr_pval_t

    cat ("Replicate: ", s, "mrtrio pval = ", pval[s],"\n") 
    cat ("Replicate: ", s, "brumpton pval  = ", pval2[s],"\n") 
    cat ("Replicate: ", s, "ivw pval  = ", pval3[s],"\n") 
    cat ("Replicate: ", s, "egger pval  = ", pval4[s],"\n") 
    cat ("Replicate: ", s, "median pval  = ", pval5[s],"\n") 
    cat ("Replicate: ", s, "mode pval  = ", pval6[s],"\n") 
    cat ("Replicate: ", s, "Brumpton-orig pval  = ", pval7[s],"\n") 
    cat ("Replicate: ", s, "ivw-filter pval  = ", pval10[s],"\n") 
    cat ("Replicate: ", s, "egger-filter pval  = ", pval11[s],"\n") 
    cat ("Replicate: ", s, "median-filter pval  = ", pval12[s],"\n") 
    cat ("Replicate: ", s, "mode-filter pval  = ", pval13[s],"\n") 
    # cat ("Replicate: ", s, "Brumpton-tdist pval  = ", pval8[s],"\n") 
    # cat ("Replicate: ", s, "IVW-tdist pval  = ", pval9[s],"\n") 

}

#Test calibration
cal_05=get_calibration(pval, 0.05)
cal2_05=get_calibration(pval2, 0.05)
cal3_05=get_calibration(pval3, 0.05)
cal4_05=get_calibration(pval4, 0.05)
cal5_05=get_calibration(pval5, 0.05)
cal6_05=get_calibration(pval6, 0.05)
cal7_05=get_calibration(pval7, 0.05)

cal10_05=get_calibration(pval10, 0.05)
cal11_05=get_calibration(pval11, 0.05)
cal12_05=get_calibration(pval12, 0.05)
cal13_05=get_calibration(pval13, 0.05)
# cal8_05=get_calibration(pval8, 0.05)
# cal9_05=get_calibration(pval9, 0.05)


cat ("## mrtrio Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal_05$n0, cal_05$n, cal_05$p,cal_05$pp,"\n")

cat ("## Brumpton Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal2_05$n0, cal2_05$n, cal2_05$p, cal2_05$pp,"\n")

cat ("## IVW Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal3_05$n0, cal3_05$n, cal3_05$p, cal3_05$pp,"\n")

cat ("## Egger Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal4_05$n0, cal4_05$n, cal4_05$p, cal4_05$pp,"\n")

cat ("## Median Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal5_05$n0, cal5_05$n, cal5_05$p, cal5_05$pp,"\n")

cat ("## Mode Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal6_05$n0, cal6_05$n, cal6_05$p, cal6_05$pp,"\n")

cat ("## Brumpton-orig Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
cat ("##",ext_model, trio_model, cal7_05$n0, cal7_05$n, cal7_05$p, cal7_05$pp,"\n")



# cat ("## IVW-filter Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
# cat ("##",ext_model, trio_model, cal10_05$n0, cal10_05$n, cal10_05$p, cal10_05$pp,"\n")

# cat ("## Egger-filter Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
# cat ("##",ext_model, trio_model, cal11_05$n0, cal11_05$n, cal11_05$p, cal11_05$pp,"\n")

# cat ("## Median-filter Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
# cat ("##",ext_model, trio_model, cal12_05$n0, cal12_05$n, cal12_05$p, cal12_05$pp,"\n")

# cat ("## Mode-filter Ext_model trio_model n0 n FPR@0.05 p-value_of_calibration\n")
# cat ("##",ext_model, trio_model, cal13_05$n0, cal13_05$n, cal13_05$p, cal13_05$pp,"\n")
