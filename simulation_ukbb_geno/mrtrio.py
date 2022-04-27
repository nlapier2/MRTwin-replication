import sys
import generate_offspring
import mr_functions
import numpy as np


def compute_quantile(child_stat, twin_stats):
	"""
	Function to compute the quantile of the child test statistic among the null
		distribution of the digital twin test statistics, e.g. the p-value.

	Args
	-----
	child_stat: MR-Trio test statistic for the true children of the trios
	twin_stats: MR-Trio test statistic for the simulated digital twins

	Returns
	-----
	quantile of child's test statistic among the twin statistics, a.k.a. the p-value
	"""
	num_less = sum([child_stat <= stat for stat in twin_stats])
	quantile = float(1 + num_less) / float(len(twin_stats) + 1)
	return quantile


def mrtrio(beta_e, se_e, beta_o, se_o, outcome_trait, par1_geno=None, par2_geno=None, child_geno=None,
			sib_genos=None, num_twins=999, mr_method='ivw'):
	"""
	Main external function for computing the MR-Trio p-value.

	Args
	-----
	par1_geno: genotypes for one set of parents for the trios or duos
	par2_geno: genotypes for the other set of parents for the trios
	child_geno: genotypes for the children of the trios or single parents
	sib_genos: genotypes for siblings in sibling mode
	outcome_trait: list containing the outcome trait values for each child
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait
	se_e:   numpy array of the standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the standard error associated with
			beta_o
	num_twins: number of digital twins to simulate
	mr_method: which MR method to compute the statistic based on
				(ivw, egger, median, or mode)

	Returns
	-----
	res: results dictionary containing MR IVR and MR Trio results
	"""
	# check for valid input arguments
	if (par1_geno is None or child_geno is None) and sib_genos is None:
		sys.exit('Error: must supply par1_geno and child_geno (plus optionally par2_geno) or sib_genos.')
	elif (par1_geno is not None or par2_geno is not None or child_geno is not None) and sib_genos is not None:
		sys.exit('Error: both parents/children and siblings provided -- ambiguous mode. Provide one or the other.')
	elif sib_genos is not None and len(sib_genos) == 1:
		sys.exit('Error: multiple siblings required in sibling mode.')
	mr_method = mr_method.lower()
	if mr_method not in ['ivw', 'egger', 'median', 'mode']:
		sys.exit('Error: acceptable mr_method options are ivw, egger, median, or mode')

	# obtain standard (non-robust) MR causal effect estimate based on the provided betas and stderrs
	beta_e, se_e, beta_o, se_o = np.array(beta_e), np.array(se_e), np.array(beta_o), np.array(se_o)
	if mr_method == 'egger':
		correlation_matrix = None
		causal_est, intercept_est = mr_functions.mr_egger(beta_e, se_e, beta_o, se_o, correlation_matrix)
		mr_est, mr_se, mr_stat, mr_pval = causal_est
	elif mr_method == 'median':
		mr_est, mr_se, mr_stat, mr_pval = mr_functions.mr_median(beta_e, se_e, beta_o, se_o)
	elif mr_method == 'mode':
		mr_est, mr_se, mr_stat, mr_pval = mr_functions.mr_mode(beta_e, se_e, beta_o, se_o)
	else:  # mr_method == 'ivw':
		mr_est, mr_se, mr_stat, mr_pval = mr_functions.mr_ivw(beta_e, se_e, beta_o, se_o)

	# simulate digital twins, compute statistics for real & synthetic children/sibs, and get permuted p-value
	twin_stats = []
	if sib_genos is None:  # trio mode or parent-child duo mode, depending on whether par2_geno is provided
		real_stat = mr_functions.compute_statistic(
			beta_e, se_e, beta_o, se_o, child_geno, outcome_trait, mr_method, mr_est)
		for i in range(num_twins):
			counts_child_geno = None
			if par2_geno is not None:  # trio mode
				twin_genos = generate_offspring.generate_offspring(par1_geno, par2_geno, 1)
			else:  # parent-child duo mode
				twin_genos, counts_child_geno = generate_offspring.generate_offspring_duo(
					par1_geno, child_geno, counts_child_geno=counts_child_geno)
			twin_stats.append(mr_functions.compute_statistic(
				beta_e, se_e, beta_o, se_o, twin_genos, outcome_trait, mr_method, mr_est))
	else:
		real_stat = mr_functions.compute_statistic_sib(
			beta_e, se_e, beta_o, se_o, sib_genos, outcome_trait, mr_method, mr_est)
		# twin_prob_dict, count_dict = None, None
		for i in range(num_twins):
			twin_genos = generate_offspring.generate_digital_sibs_shuf(sib_genos)
			twin_stats.append(mr_functions.compute_statistic_sib(
				beta_e, se_e, beta_o, se_o, twin_genos, outcome_trait, mr_method, mr_est))
	mrtwin_p_value = compute_quantile(real_stat, twin_stats)

	# print and return results
	res = {'MR Estimate': mr_est, 'MR std. error': mr_se, 'MR t-statistic': mr_stat,
			'MR p-value': mr_pval, 'MR Twin p-value': mrtwin_p_value}
	printable_res = {'MR Estimate': "{:.2e}".format(mr_est),
						'MR std. error': "{:.2e}".format(mr_se),
						'MR t-statistic': "{:.2e}".format(mr_stat),
						'MR p-value': "{:.2e}".format(mr_pval),
						'MR Twin p-value': "{:.2e}".format(mrtwin_p_value)}
	print(printable_res)
	return res
