import argparse
import numpy as np
from bed_reader import open_bed, sample_file
import statsmodels.api as sm


def parseargs():  # handle user arguments
    parser = argparse.ArgumentParser(description='Regress a phenotype on PLINK BED genotypes.')
    parser.add_argument('--bfile', required=True, help='Plink bfile base name. Required.')
    parser.add_argument('--pheno', required=True, help='Plink phenotype file. Required.')
    parser.add_argument('--snps', required=True, help='List of pruned SNPs to use.')
    parser.add_argument('--output', default='gwas_results.tsv', help='Output file name. Default: gwas_results.tsv')
    parser.add_argument('--std_geno', action='store_true', help='Use to standardize genotypes.')
    parser.add_argument('--std_pheno', action='store_true', help='Use to standardize phenotypes.')
    user_args = parser.parse_args()
    return user_args


def read_snps(snpfile):
    # read in a list of pruned snps to keep
    snp_dict, snp_list = {}, []
    with(open(snpfile, 'r')) as infile:
        for line in infile:
            snp = line.strip()
            snp_dict[snp] = True
            snp_list.append(snp)
    return snp_dict, snp_list


def find_snp_positions(bfile, snp_list):
    # find ordered place of snps in snp_list in the bimfile so we know which snps to read later
    pos = []
    snp2a1 = {}
    bimfile = bfile + '.bim'
    with(open(bimfile, 'r')) as infile:
        counter = 0
        for line in infile:
            splits = line.strip().split()
            snp, a1 = splits[1], splits[4]
            if snp in snp_list:
                pos.append(counter)
                snp2a1[snp] = a1
            counter += 1
    return pos, snp2a1


def read_fam_iids(fname):
    # read iids from a plink fam file
    iids = {}
    counter = 0
    with(open(fname, 'r')) as infile:
        for line in infile:
            iids[line.split()[0]] = counter
            counter += 1
    return iids


def read_pheno_iids(fname):
    # read iids from a plink pheno file
    iids = []
    with(open(fname, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            iids.append(line.split()[0])
    return iids


def get_reordering(fam_iids1, pheno_iids1):
    # figure out how to reorder fam people to be in same order as pheno file
    reordering = []
    for iid in pheno_iids1:
        fam_order = fam_iids1[iid]
        reordering.append(fam_order)
    return np.array(reordering)


def read_pheno(fname, std_pheno):
    all_phenos = []
    with(open(fname, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            pheno = float(line.strip().split()[2])
            all_phenos.append(pheno)
    all_phenos = np.array(all_phenos)
    if std_pheno:
        all_phenos = (all_phenos - np.mean(all_phenos)) / np.std(all_phenos)
    return all_phenos


def gwas_pass(bfile, reordering, outname, pheno, std_geno, snp_positions, snp_list, snp2a1):
    bed = open_bed(bfile + '.bed')
    with(open(outname, 'w')) as outfile:
        outfile.write('.\t.\tSNP_ID\tA1\t.\t.\t.\t.\tBETA\tSTDERR\tP-VALUE\n')
        # read in the genotypes from the desired snps
        all_genos = []
        for pos in snp_positions:
            genos = bed.read(index=np.s_[:, pos]).T[0]
            if std_geno:
                genos = (genos - np.mean(genos)) / np.std(genos)
            genos_reord = genos[reordering]
            all_genos.append(genos_reord)
        all_genos = np.array(all_genos).T

        # run regression and write results
        ols_res = sm.OLS(pheno, all_genos, missing='drop').fit()
        for i in range(len(snp_list)):
            snp_id = snp_list[i]
            a1 = snp2a1[snp_id]
            beta = ols_res.params[i]
            stderr = ols_res.bse[i]
            pval = ols_res.pvalues[i]
            beta, stderr, pval = str(beta), str(stderr), str(pval)
            outfile.write('\t'.join(['.', '.', snp_id, a1, '.', '.', '.', '.', beta, stderr, pval]) + '\n')


if __name__ == '__main__':
    args = parseargs()
    # find which lines in the bfile to read in, corresponding to SNPs in the args.snps file
    snpdict, snplist = read_snps(args.snps)
    snp_pos, snp_to_a1 = find_snp_positions(args.bfile, snpdict)

    # determine how to reorder snps so genos and phenos are in same order
    fam_iids = read_fam_iids(args.bfile + '.fam')
    pheno_iids = read_pheno_iids(args.pheno)
    reord = get_reordering(fam_iids, pheno_iids)

    # read phenotypes and run gwas
    phenos = read_pheno(args.pheno, args.std_pheno)
    gwas_pass(args.bfile, reord, args.output, phenos, args.std_geno, snp_pos, snplist, snp_to_a1)

