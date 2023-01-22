import argparse
import sys
sys.path.insert(1, "../MR-Twin")
import mrtwin
import numpy as np


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Evaluate MR-Twin & standard MR on simulation results.')
    parser.add_argument('--basename', required=True, help='Base name for input files. Required.')
    parser.add_argument('--num_twins', type=int, default=100, help='Number of digital twins to simulate. Default: 100')
    parser.add_argument('--method', default='ivw', choices=['ivw', 'egger', 'median', 'mode'],
                        help='Use ivw, egger, median, or mode statistic. Default: ivw.')
    parser.add_argument('--include_duo', action='store_true', help='Include duo mode in the analysis.')
    parser.add_argument('--include_sib', action='store_true', help='Include duo mode in the analysis.')
    # parser.add_argument('--nsims', type=int, help='number of simulations performed')
    # parser.add_argument('--include_sib', action='store_true', help='Include sibling mode in the analysis.')
    user_args = parser.parse_args()
    return user_args


def read_infiles(basename):
    betas_and_stderrs = np.loadtxt(basename + '.betahatEO')
    beta_e, se_e, beta_o, se_o = betas_and_stderrs
    #print(se_o)
    #print(se_e)
    pa1geno = np.loadtxt(basename + '.pa1geno')
    pa2geno = np.loadtxt(basename + '.pa2geno')
    child1_geno = np.loadtxt(basename + '.childgeno')
    child1_pheno = np.loadtxt(basename + '.childpheno')
    child1_e, child1_o = np.transpose(child1_pheno[0]), np.transpose(child1_pheno[1])
    return beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o


def read_sib_infiles(basename):
    betas_and_stderrs = np.loadtxt(basename + '.betahatEO')
    beta_e, se_e, beta_o, se_o = betas_and_stderrs
    #print(se_o)
    #print(se_e)
    pa1geno = np.loadtxt(basename + '.pa1geno')
    pa2geno = np.loadtxt(basename + '.pa2geno')
    child1_geno = np.loadtxt(basename + '.childgeno')
    child1_pheno = np.loadtxt(basename + '.childpheno')
    
    child2_geno = np.loadtxt(basename + '.child2geno')
    child2_pheno = np.loadtxt(basename + '.child2pheno')
    child1_e, child1_o = np.transpose(child1_pheno[0]), np.transpose(child1_pheno[1])
    child2_e, child2_o = np.transpose(child2_pheno[0]), np.transpose(child2_pheno[1])
    return beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o, child2_geno,child2_o

if __name__ == '__main__':
    main_args = parseargs()
    beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o = read_infiles(main_args.basename)
    #print(se_e)
    #print(se_o)
    #print(beta_e)
    #print(beta_o)
    print("## MR-Twin:")
    mr_res = mrtwin.mrtwin(beta_e,se_e,beta_o, se_o,child1_o, pa1geno, pa2geno, child1_geno,
                           num_twins=main_args.num_twins, mr_method=main_args.method)
    
    
    if main_args.include_duo:
        print("## MR-Duo:")
        mr_duo_res = mrtwin.mrtwin(beta_e, se_e, beta_o, se_o, child1_o, par1_geno=pa1geno, child_geno=child1_geno,
                                       num_twins=main_args.num_twins, mr_method=main_args.method)
        
    
    if main_args.include_sib:
        print("## MR-Sib:")
        beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o, child2_geno, child2_o = read_sib_infiles(main_args.basename)
        sib_genos=np.array([child1_geno,child2_geno])
        outcomes = [child1_o,child2_o]
        mr_duo_res = mrtwin.mrtwin(beta_e, se_e, beta_o, se_o, outcomes, par1_geno=None, par2_geno=None,child_geno=None,
                                       num_twins=main_args.num_twins, mr_method=main_args.method, sib_genos=sib_genos)
