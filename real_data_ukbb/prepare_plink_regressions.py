# Given a dataframe of features and a plink fam file, prepare pheno and covar files to run plink GWAS on
import argparse
import sys


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Prepare pheno and covar files for PLINK.')
    parser.add_argument('--feature_file', required=True, help='Input file with phenotypes and PCs. Required.')
    parser.add_argument('--famfile', required=True, help='PLINK fam file. Required.')
    parser.add_argument('--outname', required=True, help='Base name for pheno and covar files. Required.')
    parser.add_argument('--pheno', required=True, help='Name of phenotype to extract from feature file.')
    user_args = parser.parse_args()
    return user_args


def read_famfile(fname):
    iid_dict = {}
    with(open(fname, 'r')) as infile:
        for line in infile:
            iid = line.split(' ')[1]
            iid_dict[iid] = True
    return iid_dict


def get_pheno_covar_cols(pheno, covars, line):
    # determine which columns correspond to target phenotype and covariates
    pheno_col, iid_col = -1, -1
    covar_cols = []
    covar_list = []  # tracks order of covariates, since covars is a dict
    header = line.strip().split('\t')
    for i in range(len(header)):
        if header[i] == pheno:
            pheno_col = i
        elif header[i] in covars:
            covar_cols.append(i)
            covar_list.append(header[i])
        elif header[i] == 'EID':
            iid_col = i
    if pheno_col == -1 or iid_col == -1:
        sys.exit('Error: pheno or EID column not found.')
    return pheno_col, iid_col, covar_cols, covar_list


def write_pheno_covar(iid_dict, feature_file, outname, target_pheno, target_covars):
    with(open(feature_file, 'r')) as infile:
        with(open(outname + '.pheno', 'w')) as pout:
            with(open(outname + '.covar', 'w')) as cout:
                pheno_col, iid_col, covar_cols, covar_list = get_pheno_covar_cols(
                    target_pheno, target_covars, infile.readline())
                pout.write('FID IID pheno\n')  # pheno file header line
                cout.write('FID IID ' + ' '.join(covar_list) + '\n')  # covar file header line
                for line in infile:
                    splits = line.strip().split('\t')
                    iid, pheno = splits[iid_col], splits[pheno_col]
                    if iid not in iid_dict or pheno == '':
                        continue
                    pout.write(iid + ' ' + iid + ' ' + pheno + '\n')
                    covars = [splits[i] for i in covar_cols]
                    for i in range(len(covars)):
                        if covars[i] == 'Male':
                            covars[i] = '0'
                        if covars[i] == 'Female':
                            covars[i] = '1'
                    covars = ' '.join(covars)
                    cout.write(iid + ' ' + iid + ' ' + covars + '\n')


if __name__ == '__main__':
    args = parseargs()  # parse user-passed arguments
    # define covariates
    main_covars = {'age': True, 'sex': True}
    num_pcs = 20
    for z in range(1, num_pcs + 1):
        main_covars['PC' + str(z)] = True

    # read in IIDs from fam file, read features for those IIDs, and write to pheno and covar files
    main_iid_dict = read_famfile(args.famfile)
    write_pheno_covar(main_iid_dict, args.feature_file, args.outname, args.pheno, main_covars)
