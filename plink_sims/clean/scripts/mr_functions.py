import numpy as np
import statsmodels.api as sm
from scipy import stats, optimize
import sys
from KDE import BinDist, KDE


def mr_ivw(beta_e, se_e, beta_o, se_o):
	"""
	MR IVW implementation
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait
	se_e:   numpy array of the  standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the  standard error associated with
			beta_o
	Returns
	-----
	mr_estimate, standard error associated with mr_est and
	mr_stat(chisq) with the associated p-value
	"""
	
	model = sm.WLS(beta_o, beta_e, weights=np.divide(1, np.square(se_o)))
	res = model.fit()
	est = res.params[0]
	se = res.bse[0]/min(1,(res.scale))
	t = est/se
	p = 2*stats.t.cdf(-abs(t), df = len(beta_e) - 1)
	return est,se,t,p


def mr_mode_est(ratio_est, ratio_se, phi):
	"""
	Estimation of mode-based single casual estimate 
	for multiple genetic instruments 
	Args
	-----
	ratio_est: Observed casual effects of each variant (Beta_iv)
	ratio_se: standard error of observed casual effects of each variant
	phi: Bandwidth parameter
	Returns
	-----
	Mode-based estimate
	"""
	s = 0.9 * (min(np.std(ratio_est), stats.median_abs_deviation(ratio_est, scale='normal'))) / (
				len(ratio_est) ** (1 / 5))
	weights_s = (np.reciprocal(np.square(ratio_se))) / (np.sum(np.reciprocal(np.square(ratio_se))))
	h = s * phi
	#kde = stats.gaussian_kde(ratio_est, weights=weights_s, bw_method=h)
	#opt = optimize.minimize_scalar(lambda x: - kde(x))
	#beta_est = opt.x
	x_a, y_a = KDE(ratio_est,bw = h, weights = weights_s)
	beta_est = x_a[max(y_a) == y_a]
	#beta_est = subprocess.call(['density.R', ratio_est, h, weights_s])
	return beta_est


def mr_mode_boot(ratio_est, ratio_se, weighted, it=10000, phi=1):
	"""
	Employs the bootstrap technique to estimate the standard error of
	the MR mode-based estimate
	Args
	-----
	ratio_est: Observed casual effects of each variant (Beta_iv)
	ratio_se: standard error of observed casual effects of each variant
	weighted: weights for SNPs
	phi: Bandwidth parameter
	it: number of iterations for the bootstrap technique
	Returns
	-----
	Estimate of the standard error of MR mode-based estimate
	"""
	betas_bootstrap = [0] * it
	for i in range(it):
		beta_iv = np.random.normal(loc=ratio_est, scale=ratio_se)
		if weighted:
			betas_bootstrap[i] = mr_mode_est(beta_iv, ratio_se, phi)
		else:
			betas_bootstrap[i] = mr_mode_est(beta_iv, [1] * len(ratio_est), phi)
	return betas_bootstrap


def mr_mode(beta_e, se_e, beta_o, se_o, weighted=True, phi=1, it=10000):
	"""
	MR mode-based estimation implementation
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait 
	se_e:   numpy array of the  standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the  standard error associated with
			beta_o
	weighted: weights for SNPs
	phi: Bandwidth parameter
	it: number of iterations for bootstrap
	Returns
	-----
	MR mode-based estimate, MR mode-based standard-error, 
	corresponding z-statisitc and p-value
	"""
	beta_iv = np.divide(beta_o, beta_e)
	part_sum = np.square(beta_o) * np.divide(np.square(se_e), np.square(np.square(beta_e)))
	beta_iv_se = np.sqrt(np.divide(np.square(se_o), np.square(np.abs(beta_e))) + part_sum)
	#print(part_sum)
	#print(beta_iv_se)
	if weighted:
		beta_mode = mr_mode_est(beta_iv, beta_iv_se, phi)
	else:
		beta_mode = mr_mode_est(beta_iv, [1] * len(beta_e), phi)
	beta_mode_se = stats.median_abs_deviation(mr_mode_boot(beta_iv, beta_iv_se, weighted, it, phi), scale='normal')
	if "numpy.ndarray" in str(type(beta_mode_se)):
		beta_mode_se = float(beta_mode_se)
	z_stat = beta_mode / beta_mode_se
	return beta_mode[0], beta_mode_se, z_stat[0], (2 * stats.norm.cdf(-abs(z_stat)))[0]


def mr_egger(beta_e, se_e, beta_o, se_o, correlation_matrix=None):
	"""
	MR Egger implementation
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait 
	se_e:   numpy array of the  standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the  standard error associated with
			beta_o
	correlation_matrix: correlation matrix of the variants
	Returns
	-----
	Two one-dimensional 4-element lists:
	List one: Causal test statistics/estimate. Elements: Causal estimate, Causal estimate se,
	causal estimate test stat, p-value
	List two: Directional pleiotropy test. Elements: intercept estimate, intercept estimate se,
	intercept estimate test stat, p-value
	"""
	#print("do")
	if len(beta_e) < 3:
		sys.exit("More than 2 variants required")
	beta_ea = np.abs(beta_e)
	sign = np.copy(beta_e)
	sign[beta_e >= 0] = 1
	sign[sign < 0] = -1
	beta_oa = sign * beta_o
	if correlation_matrix is not None:
		rho = correlation_matrix * (np.outer(sign, np.transpose(sign)))
		omega = np.outer(se_o, np.transpose(se_o)) * rho
		thetas_p1 = np.transpose(np.array([np.ones(len(beta_ea)), beta_ea]))
		thetas_p2 = np.linalg.pinv(omega)
		thetas = np.dot(np.linalg.pinv(np.dot(np.transpose(thetas_p1), np.dot(thetas_p2, thetas_p1))), np.dot(np.transpose(thetas_p1), np.dot(thetas_p2, beta_oa)))
		theta_e = thetas[1]
		theta_i = thetas[0]
		rse = beta_oa - theta_i - theta_e * beta_ea
		rse_corr = float(np.sqrt(np.dot(np.transpose(rse), np.dot(thetas_p2, np.divide(rse, len(beta_ea) - 2)))))
		sigma = np.linalg.pinv(np.dot(np.transpose(thetas_p1), np.dot(thetas_p2, thetas_p1))) * max(
			np.sqrt(np.dot(np.transpose(rse), np.dot(thetas_p2, np.divide(rse, (len(beta_ea) - 2))))), 1)
		theta_e_se = np.sqrt(sigma[1][1]) / min(1., rse_corr)
		theta_i_se = np.sqrt(sigma[0][0]) / min(1., rse_corr)
		z_stat = theta_e / theta_e_se
		z_stat2 = theta_i / theta_i_se
	else:
		#print(beta_oa)
		#print(beta_ea)
		lm_w = sm.regression.linear_model.WLS(beta_oa, sm.add_constant(beta_ea), weights=np.divide(1, np.square(se_o)))
		results = lm_w.fit()
		rse = np.sqrt(results.mse_resid)
		#print(results.params[0])
		theta_e = results.params[1]
		theta_e_se = results.bse[1] / min(rse, 1)
		theta_i = results.params[0]
		theta_i_se = results.bse[0] / min(rse, 1)
		z_stat = theta_e / theta_e_se
		z_stat2 = theta_i / theta_i_se

	return [theta_e, theta_e_se, z_stat, 2 * stats.norm.cdf(-abs(z_stat))], [theta_i, theta_i_se, z_stat2, 2 * stats.norm.cdf(-abs(z_stat2))]


def mr_median_est(theta, weights, weighted=False):
	"""
	Calculation of MR median-based estimate
	Args
	-----
	theta: Outcome beta / exposure beta
	weights: weights for SNPs in MR median
	weighted: whether or not to use weighted method
	Returns
	-----
	Mr median-based estimate
	"""
	theta_s = np.sort(theta)
	if weighted:
		o = theta.argsort()
		weights_s = weights[o]
	else:
		weights_s = weights
	cuml = (np.cumsum(weights_s) - np.multiply(0.5, weights_s)) / sum(weights_s)
	k = sum([int(x) for x in (cuml < 0.5)])
	ratio = (0.5 - cuml[k - 1]) / (cuml[k] - cuml[k - 1])
	estimate = theta_s[k - 1] + (theta_s[k] - theta_s[k - 1]) * ratio
	return estimate


def mr_median_bootstrap_se(beta_e, se_e, beta_o, se_o, weights, it=10000, weighted=False):
	"""
	Bootstrap method to find MR median-based standard error
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait 
	se_e:   numpy array of the  standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the  standard error associated with
			beta_o
	weights: weights for SNPs in MR median
	weighted: whether or not to use weighted method
	it: number of iterations for bootstrap
	Returns
	-----
	MR median-based standard error
	"""
	med = [0] * it
	for i in range(it):
		be_boot = np.random.normal(loc=beta_e, scale=se_e)
		bo_boot = np.random.normal(loc=beta_o, scale=se_o)
		theta_boot = np.divide(bo_boot, be_boot)
		med[i] = mr_median_est(theta_boot, weights, weighted)
	return np.std(med)


def mr_median(beta_e, se_e, beta_o, se_o, it=10000, weighted=True):
	"""
	MR median-based implementation
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait 
	se_e:   numpy array of the  standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the  standard error associated with
			beta_o
	it: number of iterations for MR median standard error bootstrap
	weighted: Whether or not to use weighted method, if False equal weights will be applied
	Returns
	-----
	MR median-based estimate, MR median-based se, z-statistic for median method, Median method p value
	"""
	if len(beta_e) < 3:
		sys.exit("More than 2 variants required")
	theta = np.divide(beta_o, beta_e)
	estimate = 0
	if weighted:
		weights = np.square(np.divide(beta_e, se_o))
		estimate = mr_median_est(theta, weights, weighted)
		median_se = mr_median_bootstrap_se(beta_e, se_e, beta_o, se_o, weights, it, weighted=True)
	else:
		weights = [1 / len(beta_e)] * len(beta_e)
		estimate = mr_median_est(theta, weights, weighted)
		median_se = mr_median_bootstrap_se(beta_e, se_e, beta_o, se_o, weights, it, weighted=False)
	median_stat = estimate / median_se
	median_p = 2 * stats.norm.cdf(-abs(median_stat))
	return estimate, median_se, median_stat, median_p


def compute_statistic(beta_e, se_e, beta_o, se_o, genotypes, outcome_trait, mr_method, mr_est=None):
	"""
	Computes the MR-Trio negative squared error test statistic. Outcome values are predicted
		via a Two-Stage Least Squares (2SLS) style approach using the input genotypes, and
		the negative squared error against the true outcome values is computed and returned.
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait
	se_e:   numpy array of the standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the standard error associated with
			beta_o
	genotypes: matrix of genotypes with each row being a person
				and each column being a SNP
	outcome_trait: list containing the outcome trait values for each person
	mr_method: which MR method to compute the statistic based on
				(ivw, egger, median, or mode)
	mr_est: optionally pass in pre-computed mr estimate
	Returns
	-----
	neg_sq_error: the negative squared error statistic used for the MR-Trio test
	"""
	#print(beta_e[:10]) ; print(genotypes[0][:10])
	#print(sum(beta_e * genotypes[0]))
	est_exposures = np.sum(beta_e * genotypes, axis=1)
	if mr_est is None:
		if mr_method == 'egger':
			correlation_matrix = None
			causal_est, intercept_est = mr_egger(beta_e, se_e, beta_o, se_o, correlation_matrix)
			mr_est, mr_se, mr_stat, mr_pval = causal_est
		elif mr_method == 'median':
			mr_est, mr_se, mr_stat, mr_pval = mr_median(beta_e, se_e, beta_o, se_o)
		elif mr_method == 'mode':
			mr_est, mr_se, mr_stat, mr_pval = mr_mode(beta_e, se_e, beta_o, se_o)
		else:  # mr_est == 'ivw':
			mr_est, mr_se, mr_stat, mr_pval = mr_ivw(beta_e, se_e, beta_o, se_o)
		# mr_est, mr_se, mr_stat, mr_pval = mr_ivw(beta_e, se_e, beta_o, se_o)
	#est_exposures = (est_exposures - np.mean(est_exposures)) #/ np.std(est_exposures)
	pred_outcomes = est_exposures * mr_est
	###
	#outcome_trait = (outcome_trait - np.mean(outcome_trait)) #/ np.std(outcome_trait)
	#print(pred_outcomes[:10]) ; print(outcome_trait[:10])
	#pred_outcomes = pred_outcomes - np.mean(pred_outcomes)
	errs = pred_outcomes - outcome_trait
	#errs = errs**2
	#print('Min/Max/Mean/Median/Var:')
	#print(min(errs), max(errs), np.mean(errs), np.median(errs), np.var(errs))
	###print(pred_outcomes - outcome_trait)
	###
	neg_sq_error = -sum(np.square(pred_outcomes - outcome_trait))
	#outcome_trait = outcome_trait - np.mean(outcome_trait)
	#neg_sq_error = -np.median(np.square(pred_outcomes - outcome_trait))
	#neg_sq_error = stats.pearsonr(pred_outcomes, outcome_trait)[0]**2
	#print(neg_sq_error)
	return neg_sq_error


def compute_statistic_sib(beta_e, se_e, beta_o, se_o, genotypes, outcome_trait, mr_method, mr_est=None):
	"""
	Computes the MR-Trio negative squared error test statistic (above) for each set of siblings and returns the average.
	Args
	-----
	beta_e: numpy array of beta-coefficient values
			for genetic associations with the exposure trait
	se_e:   numpy array of the standard error associated with
			beta_e
	beta_o: numpy array of the beta-coefficient values
			for genetic association with the outcome
			trait
	se_o: numpy array of the standard error associated with
			beta_o
	genotypes: 3-D matrix of genotypes, with each entry being a genotype matrix for a set of sibs,
				with each row being a person and each column being a SNP
	outcome_trait: 2-D matrix, each entry containing the outcome trait values for the n-th set of sibs
	mr_method: which MR method to compute the statistic based on
				(ivw, egger, median, or mode)
	Returns
	-----
	combined_stat: the average result from running compute_statistic on each set of siblings
	"""
	stats_for_each_sib_set = []
	for i in range(len(genotypes)):
		sib_set_genotypes = genotypes[i]
		sib_set_outcome_trait = outcome_trait[i]
		stat = compute_statistic(
			beta_e, se_e, beta_o, se_o, sib_set_genotypes, sib_set_outcome_trait, mr_method, mr_est=mr_est)
		stats_for_each_sib_set.append(stat)
	combined_stat = np.mean(stats_for_each_sib_set)
	return combined_stat


def mr_all_methods(beta_e, se_e, beta_o, se_o, weighted_p=True, it_p=10000, correlation_matrix=None, phi=1):
	ivw = mr_ivw(beta_e, se_e, beta_o, se_o)
	median = mr_median(beta_e, se_e, beta_o, se_o, it=it_p, weighted=weighted_p)
	egger = mr_egger(beta_e, se_e, beta_o, se_o, correlation_matrix)
	mode = mr_mode(beta_e, se_e, beta_o, se_o, weighted=weighted_p, phi=phi, it=it_p)
	return ivw, median, egger, mode
