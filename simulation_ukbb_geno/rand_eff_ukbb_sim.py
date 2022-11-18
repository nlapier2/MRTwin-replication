# Randomly select causal SNPs from Plink bfile, sample effect sizes,
#     generate phenotypes, and run GWAS. Write output for MR methods to run on.
import argparse
import glob
import subprocess
import sys
import numpy as np
import pandas as pd
from pandas_plink import read_plink
import statsmodels.formula.api as smf
from generate_offspring import generate_offspring


# declare global-scope output file names. later set to extend --output base name.
e_causal_fname, o_causal_fname, all_causal_fname = '', '', ''
ext_iid_fname, parent_iid_fname = '', ''
ext_plink_fname, trio_plink_fname, combined_plink_fname = '', '', ''
pheno_df_fname = ''
gwas_out_fname = ''
pa1geno_fname, pa2geno_fname = '', ''
child1geno_fname, child1pheno_fname = '', ''
child2geno_fname, child2pheno_fname = '', ''
betahat_eo_fname = ''


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Random effects GWAS simulator for UK Biobank genotypes.')
    parser.add_argument('--bfile', required=True, help='Plink bfile base name. Required.')
    parser.add_argument('--evecs', required=True, help='Eigenvectors file w/ one SNP per line. Required.')
    parser.add_argument('--h2ge', type=float, default=0.25, help='Heritability of genetic effects on exposure.')
    parser.add_argument('--h2go', type=float, default=0.25, help='Heritability of genetic effects on outcome.')
    parser.add_argument('--h2ce', type=float, default=0.0, help='Heritability of confounder effect on exposure.')
    parser.add_argument('--h2co', type=float, default=0.25, help='Heritability of confounder effect on outcome.')
    parser.add_argument('--h2eo', type=float, default=0.25, help='Heritability of exposure effect on outcome.')
    parser.add_argument('--no_gwas', action='store_true', help='Use to not run GWAS on the simulated phenotypes.')
    parser.add_argument('--num_causal', type=int, default=1000, help='Number of SNPs to make causal.')
    parser.add_argument('--num_pcs_confound', type=int, default=10, help='Number of PCs to affect confounder.')
    parser.add_argument('--num_pcs_regress', type=int, default=5, help='Number of PCs to regress out.')
    parser.add_argument('--output', default='sim_results', help='Prefix for output files.')
    parser.add_argument('--trio_pct', type=float, default=0.1, help='Percent of individuals to use for trio dataset.')
    user_args = parser.parse_args()
    return user_args


# compute projections into PC space
def make_pc_df(bfile, evecs_file):
    # read eigenvectors and genotype data
    evecs = np.loadtxt(evecs_file)
    bedfile, famfile = bfile + '.bed', bfile + '.fam'
    (bim, fam, bed) = read_plink(bedfile, verbose=False)
    # project into PC space and store in dataframe
    proj = np.array(np.matmul(bed.transpose(), evecs))
    cols = ['PC' + str(i) for i in range(1, num_pcs + 1)]
    pc_df = pd.DataFrame(proj, columns=cols)
    pc_df['EID'] = fam['iid']
    return pc_df


def create_confounder_var(pc_df, num_pcs_confound):  # confounder = sum of pcs
    conf = pc_df['PC1'] * 0.0  # initialize to zeros
    for i in range(1, num_pcs_confound + 1):
        pc = pc_df['PC'+str(i)]
        conf += pc
    if num_pcs_confound > 0:
        pc_df['conf'] = (conf - np.mean(conf)) / np.std(conf)  # standardize & add to dataframe
    else:
        pc_df['conf'] = conf
    return pc_df


def sample_causal_snps(bfile, h2ge, h2go, num_causal):
    # extract snps from bfile
    snplist = []
    with(open(bfile + '.bim', 'r')) as infile:
        for line in infile:
            snplist.append(line.split('\t')[1])
    np.random.shuffle(snplist)  # in-place random shuffle

    # now we can just pick in "order" since snplist has been shuffled
    e_causal_snps = snplist[:num_causal]
    o_causal_snps = snplist[num_causal:2*num_causal]
    e_causal = {snp: np.random.normal(0, np.sqrt(h2ge / num_causal)) for snp in e_causal_snps}
    o_causal = {snp: np.random.normal(0, np.sqrt(h2go / num_causal)) for snp in o_causal_snps}
    return snplist, e_causal, o_causal


def write_causal(e_causal, o_causal):
    with(open(e_causal_fname, 'w')) as outfile:
        for snp in e_causal:
            outfile.write(snp + '\n')
    with(open(o_causal_fname, 'w')) as outfile:
        for snp in o_causal:
            outfile.write(snp + '\n')
    with(open(all_causal_fname, 'w')) as outfile:
        subprocess.Popen(['cat', e_causal_fname, o_causal_fname], stdout=outfile).wait()


# Randomly partition iids into external and parents for trios and split plink files
def split_ids_and_plink(bfile, trio_pct):
    # Collect all iids from fam file
    iid_list = []
    with(open(bfile + '.fam', 'r')) as infile:
        for line in infile:
            iid = line.split()[0]
            iid_list.append(iid)

    # now make the random picks
    np.random.shuffle(iid_list)
    num_trio_parents = int(trio_pct * len(iid_list))
    if num_trio_parents % 2 == 1:
        num_trio_parents += 1  # ensure even number of parents for trios
    parent_iids = iid_list[:num_trio_parents]
    ext_iids = iid_list[num_trio_parents:]

    # write picked trio and external iids to files, then use plink to subset original file using the iids & causal snps
    with(open(ext_iid_fname, 'w')) as outfile:
        for iid in ext_iids:
            outfile.write(iid + ' ' + iid + '\n')
    with(open(parent_iid_fname, 'w')) as outfile:
        for iid in parent_iids:
            outfile.write(iid + ' ' + iid + '\n')
    subprocess.Popen(['plink', '--bfile', bfile, '--keep', ext_iid_fname,
                      '--make-bed', '--memory', '4096', '--out', '1'+ext_plink_fname]).wait()
    subprocess.Popen(['plink', '--bfile', bfile, '--keep', parent_iid_fname,
                      '--recode', '12', '--memory', '4096', '--out', '1'+trio_plink_fname]).wait()
    return ext_iids, parent_iids


def get_genotypes(splits):  # combine plink 12 alleles into 0/1/2 genotypes
    genos = []
    for i in range(int(len(splits) / 2)):
        genotype = int(splits[2*i]) + int(splits[2*i+1]) - 2
        genos.append(float(genotype))
    return genos


def get_alleles(genos):  # split genotypes into plink12 alleles
    alleles = ''
    for g in genos:
        if g == 0:
            alleles += ' 1 1'
        elif g == 1:
            alleles += ' 2 1'
        else:
            alleles += ' 2 2'
    return alleles


def sim_children():  # simulate children for trios, add them to plink PED file, then create tfile
    iid_list = []
    pa1genos, pa2genos = [], []
    counter = 0
    # read ped file to get parent genos
    with(open(trio_plink_fname + '.ped', 'r')) as pedfile:
        for line in pedfile:
            counter = (counter + 1) % 2  # pick whether parent 1 or 2
            splits = line.strip().split()
            iid_list.append(int(splits[0]))
            genos = get_genotypes(splits[6:])
            if counter == 0:
                pa1genos.append(genos)
            else:
                pa2genos.append(genos)

    # shuffle parent genos to set up random mating, then generate children & their unique iids
    pa1genos, pa2genos = np.array(pa1genos), np.array(pa2genos)
    np.random.shuffle(pa1genos)
    np.random.shuffle(pa2genos)
    child1_geno = generate_offspring(pa1genos, pa2genos, 1)
    child2_geno = generate_offspring(pa1genos, pa2genos, 1)

    # append child genos to ped file
    # max_pa_iid = max(iid_list) + 1  # start counting child iids after max parental iid, to avoid duplicates
    child1_iid_list = [str(999999999 + i) for i in range(len(child1_geno))]
    child2_iid_list = [str(999999999 + i + len(child1_geno)) for i in range(len(child2_geno))]
    with(open(trio_plink_fname + '.ped', 'a')) as pedfile:
        static_fields = ' 0 0 1 -9'
        for i in range(len(child1_geno)):
            plink12_alleles_1 = get_alleles(child1_geno[i])
            iid_fid_1 = child1_iid_list[i] + ' ' + child1_iid_list[i]
            pedfile.write(iid_fid_1 + static_fields + plink12_alleles_1 + '\n')
            plink12_alleles_2 = get_alleles(child2_geno[i])
            iid_fid_2 = child2_iid_list[i] + ' ' + child2_iid_list[i]
            pedfile.write(iid_fid_2 + static_fields + plink12_alleles_2 + '\n')
    return child1_iid_list, child2_iid_list, pa1genos, pa2genos, child1_geno, child2_geno


def generate_genotypic_effect(bfile, causal_snps):  # generate genotypic effects on phenotypes
    genotypic_effect = []
    (bim, fam, bed) = read_plink(bfile + ".bed", verbose=False)
    for i in range(len(bed)):
        snp_id = bim[bim['i']==i]['snp'].to_string(index=False)
        if snp_id not in causal_snps:
            continue
        snp_effect = causal_snps[snp_id]
        genos = np.array(bed[i])
        genos = (genos - np.mean(genos)) / np.std(genos)  # standardize
        if len(genotypic_effect) == 0:
            genotypic_effect = snp_effect * genos
        else:
            genotypic_effect += snp_effect * genos
    return genotypic_effect


def generate_phenotypes(causal_snps, h2g, h2c, confounder, h2eo=None, ext_exposure=None, trio_exposure=None):
    ext_geno_effect = generate_genotypic_effect(ext_plink_fname, causal_snps)
    trio_geno_effect = generate_genotypic_effect(trio_plink_fname, causal_snps)
    confounding_effect = np.sqrt(h2c) * confounder
    if ext_exposure is None:  # exposure not provided, so this is the exposure trait
        ext_noise = np.random.normal(0, np.sqrt(1.0 - h2g - h2c), len(ext_geno_effect))
        trio_noise = np.random.normal(0, np.sqrt(1.0 - h2g), len(trio_geno_effect))
        ext_pheno = ext_geno_effect + confounding_effect + ext_noise
        trio_pheno = trio_geno_effect + trio_noise
    else:  # exposure provided, so this is the outcome trait
        ext_exposure_effect = np.sqrt(h2eo) * ext_exposure
        trio_exposure_effect = np.sqrt(h2eo) * trio_exposure
        ext_noise = np.random.normal(0, np.sqrt(1.0 - h2g - h2c - h2eo), len(ext_geno_effect))
        trio_noise = np.random.normal(0, np.sqrt(1.0 - h2g - h2eo), len(trio_geno_effect))
        ext_pheno = ext_geno_effect + confounding_effect + ext_exposure_effect + ext_noise
        trio_pheno = trio_geno_effect + trio_exposure_effect + trio_noise
    return ext_pheno, trio_pheno


def write_child_pheno(child1_iid_list, child2_iid_list, pheno_df):
    child1_e = [float(pc_df[pc_df['EID'] == id]['exposure']) for id in child1_iid_list]
    child1_o = [float(pc_df[pc_df['EID'] == id]['outcome']) for id in child1_iid_list]
    child2_e = [float(pc_df[pc_df['EID'] == id]['exposure']) for id in child2_iid_list]
    child2_o = [float(pc_df[pc_df['EID'] == id]['outcome']) for id in child2_iid_list]
    child1_phenos = np.array([child1_e, child1_o])
    np.savetxt(child1pheno_fname, child1_phenos)
    child2_phenos = np.array([child2_e, child2_o])
    np.savetxt(child2pheno_fname, child2_phenos)


def merge_phenos_into_df(pc_df, exposures, outcomes):
    pc_df = pc_df.sort_values(by=['EID'])
    pc_df['exposure'] = exposures
    pc_df['outcome'] = outcomes
    pc_df.to_csv(pheno_df_fname, sep='\t')
    return pc_df


def gwas_pass(bfile, pheno_df, ext_iids, num_pcs_regress):
    sig_snps, sig_snp_indices = [], []
    pc_str = ' + '.join(['PC' + str(i) for i in range(1, 1 + num_pcs_regress)])
    # read in data and extract external portion
    pheno_df = pheno_df[pheno_df['EID'].isin(ext_iids)]
    (bim, fam, bed) = read_plink(bfile + ".bed", verbose=False)
    bed = bed[:,fam[fam['iid'].isin(ext_iids)]['i']]
    
    # go through SNPs and perform linear regression including PCs on each
    counter = -1
    for i in range(len(bed)):
        counter += 1
        snp_id = bim[bim['i']==i]['snp'].to_string(index=False)
        pheno_df['genos'] = bed[i]
        if pc_str != '':
            ols_res_e = smf.ols(formula='exposure ~ genos + ' + pc_str, data=pheno_df).fit()
            ols_res_o = smf.ols(formula='outcome ~ genos + ' + pc_str, data=pheno_df).fit()
        else:  # no principal components used
            ols_res_e = smf.ols(formula='exposure ~ genos', data=pheno_df).fit()
            ols_res_o = smf.ols(formula='outcome ~ genos', data=pheno_df).fit()
        pval_e = ols_res_e.pvalues['genos']
        if pval_e < 0.001565402:  # F>10 filter; other options are 0.05 or 5*10**-8
            sig_snp_indices.append(counter)
            beta_e = ols_res_e.params['genos']
            se_e = ols_res_e.bse['genos']
            beta_o = ols_res_o.params['genos']
            se_o = ols_res_o.bse['genos']
            sig_snps.append([snp_id, beta_e, se_e, beta_o, se_o, pval_e])
    return sig_snps, sig_snp_indices


def write_gwas(sig_snps, sig_snp_indices, pa1genos, pa2genos, child1_geno, child2_geno):
    # write human-readable GWAS results
    with(open(gwas_out_fname, 'w')) as outfile:
        outfile.write('snp_id\tbeta_e\tse_e\tbeta_o\tse_o\tpval_e\n')
        for i in range(len(sig_snps)):
            outstr = '\t'.join([str(j) for j in sig_snps[i]])
            outfile.write(outstr + '\n')

    # Write beta_hats and stderr_hats for MR-Trio
    sig_snps = np.array(sig_snps)
    betahat_e = sig_snps[:, 1].astype(np.float64)
    se_e = sig_snps[:, 2].astype(np.float64)
    betahat_o = sig_snps[:, 3].astype(np.float64)
    se_o = sig_snps[:, 4].astype(np.float64)
    betahat_eo = np.array([betahat_e, se_e, betahat_o, se_o])
    np.savetxt(betahat_eo_fname, betahat_eo)

    # select significant SNPs out of trio genotypes, then write for MR-Trio (avoids weak instrument bias)
    pa1genos = pa1genos[:, sig_snp_indices]
    pa2genos = pa2genos[:, sig_snp_indices]
    child1_geno = child1_geno[:, sig_snp_indices]
    child2_geno = child2_geno[:, sig_snp_indices]
    np.savetxt(pa1geno_fname, pa1genos)
    np.savetxt(pa2geno_fname, pa2genos)
    np.savetxt(child1geno_fname, child1_geno)
    np.savetxt(child2geno_fname, child2_geno)


if __name__ == '__m__':
    args = parseargs()
    if args.bfile.endswith('.'):
        args.bfile = args.bfile[:-1]
    if any([not (0.0 <= i <= 1.0) for i in [args.h2ge, args.h2go, args.h2ce, args.h2co, args.h2eo]]):
        sys.exit('Error: all heritabilities must be between 0 and 1.')
    if args.h2ge + args.h2ce > 1.0:
        sys.exit('Error: combined heritability on exposure is greater than 1 (h2ge + h2ce > 1)')
    if args.h2go + args.h2co + args.h2eo > 1.0:
        sys.exit('Error: combined heritability on outcome is greater than 1 (h2go + h2co + h2eo > 1)')
    pre = args.output

    # set output file names
    e_causal_fname, o_causal_fname, all_causal_fname = pre + '.causal.e', pre + '.causal.o', pre + '.causal.all'
    ext_iid_fname, parent_iid_fname = pre + '.ext.iids', pre + '.parent.iids'
    ext_plink_fname, trio_plink_fname = pre + '.ext.plink', pre + '.trio.plink'
    combined_plink_fname = pre + '.comb.plink'
    pheno_df_fname = pre + '.pheno'
    gwas_out_fname = pre + '.gwas.tsv'
    pa1geno_fname, pa2geno_fname = pre + '.pa1geno', pre + '.pa2geno'
    child1geno_fname, child1pheno_fname = pre + '.child1geno', pre + '.child1pheno'
    child2geno_fname, child2pheno_fname = pre + '.child2geno', pre + '.child2pheno'
    betahat_eo_fname = pre + '.betahatEO'

    # sample causal snps, sample external & trio individuals, and sim trio children
    m_snplist, m_e_causal, m_o_causal = sample_causal_snps(args.bfile, args.h2ge, args.h2go, args.num_causal)
    write_causal(m_e_causal, m_o_causal)
    m_ext_iids, m_parent_iids = split_ids_and_plink(args.bfile, args.trio_pct)
    m_ext_iids, m_parent_iids = [int(i) for i in m_ext_iids], [int(i) for i in m_parent_iids]
    m_child1_iid_list, m_child2_iid_list, \
        m_pa1genos, m_pa2genos, m_child1_geno, m_child2_geno = sim_children()

    # merge trio & external into single bfile, then compute PC scores and generate confounder
    # write temporary 0 phenos so we have a tfam file
    subprocess.Popen(['plink', '--bfile', '1'+ext_plink_fname, '--merge', '1'+trio_plink_fname,
                          '--make-bed', '--memory', '4096', '--out', '1'+combined_plink_fname]).wait()
    m_pc_df = make_pc_df('1'+combined_plink_fname, args.evecs)
    m_pc_df = create_confounder_var(m_pc_df, args.num_pcs_confound)
    subprocess.Popen(['plink', '--bfile', '1'+combined_plink_fname, '--extract', all_causal_fname,
                      '--make-bed', '--memory', '4096', '--out', combined_plink_fname]).wait()

    # run phenotype simulation and store the results in a DataFrame with the principal components
    m_confounder = np.array(main_pc_df['conf'])
    m_exposure = generate_phenotypes(m_e_causal, args.h2ge, args.h2ce, m_confounder)
    m_outcome = generate_phenotypes(m_o_causal, args.h2go, args.h2co, m_confounder,
                                            h2eo=args.h2eo, exposure=m_exposure)
    m_pheno_df = merge_phenos_into_df(m_pc_df, m_exposure, m_outcome)
    m_ext_pheno_df = m_pheno_df[m_pheno_df['EID'].isin(m_ext_iids)]
    write_child_pheno(m_child1_iid_list, m_child2_iid_list, m_pheno_df)

    # run GWAS on the simulated phenotypes, unless the user has requested otherwise
    if not args.no_gwas:
        print('\nRunning GWAS...')
        m_sig_snps, m_sig_snp_indices = gwas_pass(ext_plink_fname, m_pheno_df, m_ext_iids, args.num_pcs_regress)
        write_gwas(m_sig_snps, m_sig_snp_indices,
                   m_pa1genos, m_pa2genos, m_child1_geno, m_child2_geno)

    # clean up files that we no longer need
    rm_cmd = ['rm'] + glob.glob(pre + '*plink*') + \
             [e_causal_fname, o_causal_fname, all_causal_fname, ext_iid_fname, parent_iid_fname, pheno_df_fname]
    subprocess.Popen(rm_cmd).wait()
