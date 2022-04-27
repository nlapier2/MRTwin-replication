import argparse
import mrtrio
import numpy as np
import os
from scipy.stats import norm


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Evaluate MR-Twin & standard MR on simulation results.')
    parser.add_argument('--dir_name', required=True, help='Directory with simulation files. Required.')
    parser.add_argument('--num_sims', type=int, default=1000, help='Number of simulations to compute results for.')
    parser.add_argument('--start_sim', type=int, default=0, help='Sim number to start on (default: 0).')
    parser.add_argument('--num_twins', type=int, default=100, help='Number of digital twins to simulate. Default: 100')
    parser.add_argument('--method', default='ivw', choices=['ivw', 'egger', 'median', 'mode'],
                        help='Use ivw, egger, median, or mode statistic. Default: ivw.')
    parser.add_argument('--include_duo', action='store_true', help='Include duo mode in the analysis.')
    parser.add_argument('--include_sib', action='store_true', help='Include sibling mode in the analysis.')
    user_args = parser.parse_args()
    return user_args


def read_sim_files(dir_name, sim_num):
    sim = str(sim_num)
    betas_and_stderrs = np.loadtxt(dir_name + sim + '.betahatEO')
    beta_e, se_e, beta_o, se_o = betas_and_stderrs
    pa1geno = np.loadtxt(dir_name + sim + '.pa1geno')
    pa2geno = np.loadtxt(dir_name + sim + '.pa2geno')
    child1_geno = np.loadtxt(dir_name + sim + '.child1geno')
    child1_pheno = np.loadtxt(dir_name + sim + '.child1pheno')
    child1_e, child1_o = np.transpose(child1_pheno[0]), np.transpose(child1_pheno[1])
    child2_geno = np.loadtxt(dir_name + sim + '.child2geno')
    child2_pheno = np.loadtxt(dir_name + sim + '.child2pheno')
    child2_e, child2_o = np.transpose(child2_pheno[0]), np.transpose(child2_pheno[1])
    return beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o, child2_geno, child2_o


def get_sim_pval(args, sim_num):
    beta_e, se_e, beta_o, se_o, pa1geno, pa2geno, child1_geno, child1_o, child2_geno, child2_o \
        = read_sim_files(args.dir_name, sim_num)
    mr_res = mrtrio.mrtrio(beta_e, se_e, beta_o, se_o, child1_o, par1_geno=pa1geno, par2_geno=pa2geno,
                           child_geno=child1_geno, num_twins=args.num_twins, mr_method=args.method)
    mr_pval = mr_res['MR p-value']
    mr_trio_pval = mr_res['MR Twin p-value']
    if args.include_duo:
        mr_duo_res = mrtrio.mrtrio(beta_e, se_e, beta_o, se_o, child1_o, par1_geno=pa1geno, child_geno=child1_geno,
                                       num_twins=args.num_twins, mr_method=args.method)
        mr_duo_pval = mr_duo_res['MR Twin p-value']
    else:
        mr_duo_pval = 1.0
    if args.include_sib:
        mr_sib_res = mrtrio.mrtrio(beta_e, se_e, beta_o, se_o, np.array([child1_o, child2_o]),
                                       sib_genos=np.array([child1_geno, child2_geno]),
                                       num_twins=args.num_twins, mr_method=args.method)
        mr_sib_pval = mr_sib_res['MR Twin p-value']
    else:
        mr_sib_pval = 1.0
    return mr_pval, mr_trio_pval, mr_duo_pval, mr_sib_pval


def get_fpr_for_setting(args):
    all_mr_pvals, all_trio_pvals, all_duo_pvals, all_sib_pvals = [], [], [], []
    completed_sims = 0
    for i in range(args.num_sims):
        this_sim_num = i + args.start_sim
        if not os.path.exists(args.dir_name + str(this_sim_num) + '.betahatEO'):
            continue
        completed_sims += 1
        mr_pval, mr_trio_pval, mr_duo_pval, mr_sib_pval = get_sim_pval(args, this_sim_num)
        all_mr_pvals.append(mr_pval)
        all_trio_pvals.append(mr_trio_pval)
        all_duo_pvals.append(mr_duo_pval)
        all_sib_pvals.append(mr_sib_pval)
    mr_fp = len([i for i in all_mr_pvals if i < 0.05])
    mr_fpr = float(mr_fp) / float(completed_sims)
    trio_fp = len([i for i in all_trio_pvals if i < 0.05])
    trio_fpr = float(trio_fp) / float(completed_sims)
    duo_fp = len([i for i in all_duo_pvals if i < 0.05])
    duo_fpr = float(duo_fp) / float(completed_sims)
    sib_fp = len([i for i in all_sib_pvals if i < 0.05])
    sib_fpr = float(sib_fp) / float(completed_sims)
    alpha = 0.05
    se = np.sqrt((alpha * (1 - alpha)) / completed_sims)  # num_sims)
    p_cal = (2 * norm.cdf(-1 * np.abs((trio_fpr - alpha) / se)))
    return mr_fpr, trio_fpr, duo_fpr, sib_fpr, p_cal, completed_sims


if __name__ == '__main__':
    main_args = parseargs()
    mr_fpr, mr_trio_fpr, mr_duo_fpr, mr_sib_fpr, p_cal, completed_sims = get_fpr_for_setting(main_args)
    print('MR Trio FPR/Power\tMR Duo FPR/Power\tMR Sib FPR/Power\tMR FPR/Power\tp-calibration')
    print(str(mr_trio_fpr) + '\t' + str(mr_duo_fpr) + '\t' + str(mr_sib_fpr) + '\t' + str(mr_fpr) + '\t' + str(p_cal))
    print('Based on ' + str(completed_sims) + ' completed simulations.')
