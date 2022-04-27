import argparse
import mrtrio
import numpy as np


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Evaluate MR-Twin & standard MR on simulation results.')
    parser.add_argument('--basename', required=True, help='Base name for input files. Required.')
    parser.add_argument('--num_twins', type=int, default=100, help='Number of digital twins to simulate. Default: 100')
    parser.add_argument('--method', default='ivw', choices=['ivw', 'egger', 'median', 'mode'],
                        help='Use ivw, egger, median, or mode statistic. Default: ivw.')
    parser.add_argument('--include_duo', action='store_true', help='Include duo mode in the analysis.')
    # parser.add_argument('--include_sib', action='store_true', help='Include sibling mode in the analysis.')
    user_args = parser.parse_args()
    return user_args


def read_infiles(basename):
    betas_and_stderrs = np.loadtxt(basename + '.betahatEO')
    beta_e, se_e, beta_o, se_o = betas_and_stderrs
    pa1geno = np.loadtxt(basename + '.pa1geno')
    pa2geno = np.loadtxt(basename + '.pa2geno')
    child1_geno = np.loadtxt(basename + '.childgeno')
    child1_pheno = np.loadtxt(basename + '.childpheno')
    child1_e, child1_o = np.transpose(child1_pheno[0]), np.transpose(child1_pheno[1])
    return beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o


if __name__ == '__main__':
    main_args = parseargs()
    beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o = read_infiles(main_args.basename)
    #print(child1_o)
    mr_res = mrtrio.mrtrio(beta_e,se_e,beta_o, se_o,child1_o, pa1geno, pa2geno, child1_geno,
                           num_twins=main_args.num_twins, mr_method=main_args.method)
    if main_args.include_duo:
        mr_duo_res = mrtrio.mrtrio(beta_e, se_e, beta_o, se_o, child1_o, par1_geno=pa1geno, child_geno=child1_geno,
                                       num_twins=main_args.num_twins, mr_method=main_args.method)
