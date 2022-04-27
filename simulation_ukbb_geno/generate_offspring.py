import numpy as np


# used in digital sibling simulation
# the nine entries correspond to the posterior probabilities of the parental genotypes being
#   '2/2, 2/1, 2/0, 1/2, 1/1, 1/0, 0/2, 0/1, or 0/0', given an observed 0, 1, or 2 genotype.
zero_posteriors = np.array([0, 0, 0, 0, 1 / 9.0, 2 / 9.0, 0, 2 / 9.0, 4 / 9.0])
one_posteriors = np.array([0, 1 / 9.0, 2 / 9.0, 1 / 9.0, 1 / 9.0, 1 / 9.0, 2 / 9.0, 1 / 9.0, 0])
two_posteriors = np.array([4 / 9.0, 2 / 9.0, 0, 2 / 9.0, 1 / 9.0, 0, 0, 0, 0])
# here we have the probability of generating a 0, 1, or 2 child for each of the 9 parent genotype combos
probs_from_parents = np.array([[1.0, 0.0, 0.0],
                               [0.5, 0.5, 0.0],
                               [0.0, 1.0, 0.0],
                               [0.5, 0.5, 0.0],
                               [0.25, 0.5, 0.25],
                               [0.0, 0.5, 0.5],
                               [0.0, 1.0, 0.0],
                               [0.0, 0.5, 0.5],
                               [0.0, 0.0, 1.0]])


def generate_offspring(father_geno, mother_geno, num_offspring, seed=None, write=False, return_parents = False, norm=False):
	"""
	Args
	-----
	father_geno: numpy array containing fathers' genotype
	mother_geno: numpy array containing mothers' genotype
	here we consider equal number of father and mother
	num_offspring: number of off-spring for each parents
	Returns
	-----
	child numpy array

	Algorithm design:
		1. generate masking sequence from bernoulli distribution (serves as a way to choose which haplotype of the parents)
		2. coin to gene: 2 -- (0/1, 1), 0 -- (0/1, 0), 1 -- (0, 0), 1 -- (1, 1)
	"""
	N = father_geno.shape[0]
	M = father_geno.shape[1]
	if seed:
		np.random.seed(seed)
	parents_geno = (father_geno,mother_geno)
	new_parents_geno = [np.zeros((0,M),float),np.zeros((0,M),float)]
	Cs = np.zeros((N,M),float)
	# for i in range(N):
	# for n in range(num_offspring):
	# C = np.zeros(M)
	for p in range(2):
		Pas = parents_geno[p].copy()
		Pa_flips = np.random.binomial(1,0.5,size=(N,M))
		Pa_1 = Pas.copy()
		Pa_1[Pas==2]=1
		Pa_1[Pas<2]=0 # mask all the other entries except 2
		Pas[Pas==2] = 0 # mask the 2 entries
		# print(f'parent geno is {parents_geno[p][i]}')
		# print(f'flips is {Pa_flips}')
		C1 = np.multiply(Pas,Pa_flips) + Pa_1
		# print(C1)
		Cs += C1
		if return_parents:
			new_parents_geno[p] = np.append(new_parents_geno[p], parents_geno[p], axis=0)
		# print(f'Child geno is {C}')
		# Cs = np.append(Cs, C.reshape(1,-1), axis=0)
	# print('done genotyping')
	if norm:
		Cs = (Cs-np.mean(Cs,axis=0))/np.std(Cs,axis=0)	
	if return_parents:
		return Cs,new_parents_geno
	else:
		return Cs


def generate_offspring_haplo(xm1, xm2, xf1, xf2, num_offspring=1, seed=None, norm=False):
	N = xm1.shape[0]
	M = xm1.shape[1]
	if seed:
		np.random.seed(seed)
	tm = np.random.binomial(1,0.5,size=(N,M))
	tf = np.random.binomial(1,0.5,size=(N,M))
	xf = np.multiply(xf1,(1-tf)) + np.multiply(xf2,tf)
	xm = np.multiply(xm1,(1-tm)) + np.multiply(xm2,tm)
	
	Cs = xm+xf
	if norm:
		Cs = (Cs-np.mean(Cs,axis=0))/np.std(Cs,axis=0)
	return Cs


def gen_duo_het_draws(parent_geno, child_geno):
	# Helper method to count number of heterozygous parent genotypes for each child genotype
	# This speeds up duo mode because we can do all np random choices at once instead of one-by-one
	counts_child_geno = {}
	for i in range(len(parent_geno)):
		for j in range(len(parent_geno[i])):
			if parent_geno[i][j] == 1:
				key = child_geno[i][j]
				if key in counts_child_geno:
					counts_child_geno[key] += 1
				else:
					counts_child_geno[key] = 1
	return counts_child_geno


def generate_offspring_duo(parent_geno, child_geno, counts_child_geno=None, seed=None, norm=False):
	"""
	Args
	-----
	parent_geno: numpy array containing parents' genotypes
	child_geno: numpy array containing children's genotypes
	seed: optionally set numpy's random seed
	norm: use to normalize the returned genotype matrix of the simulated offspring
	Returns
	-----
	child numpy array

	We fix the unknown parent's ("par2") inferred allele passed on to the child, and the digital twins randomly inherit
		one of the known parent's ("par1") alleles plus the fixed allele. Thus, if par1 is a homozygote, then both twin
		alleles are consequently fixed, so the twin will have the same genotype as the child. If par1 is a heterozygote,
		then twins will randomly inherit a 0 or 1 from them. Additionally, if the child is a homozygote, then the other
		allele (from par2) is fixed. If par1 and child are both heterozygotes, then the par1 allele is random, and the
		par2 allele is ambiguous, so the twin's genotype will be random. This can be tidily summarized as follows:

	if par1 != 1 --> twin=child
	else twin = bern(0.5) + bern(child/2)

	Algorithm:
		- Initialize the twin matrix as a copy of the child matrix
		- For any entry [i][j] of the parent matrix where par1[i][j]=1,
			set the corresponding entry in the twin matrix to bern(0.5) + bern(child[i][j]/2)
	"""
	if seed:
		np.random.seed(seed)
	if counts_child_geno is None:
		counts_child_geno = gen_duo_het_draws(parent_geno, child_geno)
	het_draws, so_far = {}, {}
	for geno in counts_child_geno:
		count = counts_child_geno[geno]
		het_draws[geno] = np.random.binomial(1, 0.5, size=count) + np.random.binomial(1, geno / 2.0, size=count)
		so_far[geno] = 0

	twin_geno = np.copy(child_geno)
	for i in range(len(parent_geno)):
		for j in range(len(parent_geno[i])):
			if parent_geno[i][j] == 1:
				# twin_geno[i][j] = np.random.binomial(1, 0.5) + np.random.binomial(1, child_geno[i][j] / 2.0)
				geno = child_geno[i][j]
				index = so_far[geno]
				twin_geno[i][j] = het_draws[geno][index]
				so_far[geno] += 1
	if norm:
		twin_geno = (twin_geno - np.mean(twin_geno, axis=0)) / np.std(twin_geno, axis=0)
	return twin_geno, counts_child_geno


def parents_given_obs(obs_genos):
	total_parent_posteriors = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])
	for geno in obs_genos:  # multiply parental posteriors over the observed sibs (assumed to be independent draws)
		if geno == 0:
			total_parent_posteriors = total_parent_posteriors * zero_posteriors
		elif geno == 1:
			total_parent_posteriors = total_parent_posteriors * one_posteriors
		else:
			total_parent_posteriors = total_parent_posteriors * two_posteriors
	total_parent_posteriors = total_parent_posteriors / sum(total_parent_posteriors)
	return total_parent_posteriors


def compute_sib_geno_prob(obs_genos):
	"""
	Args
	-----
	obs_genos: observed sibling genotypes for a particular SNP

	Returns
	-----
	geno_probs: the probability of digital sibs having 0/1/2 genotype based on the observed sibling genos

	Here is the general algorithm for computing the probabilities of [0,1,2] digital sibs for any number of siblings:
	- Initialize a numpy array mapping each of the 9 possible parental genotype pairs to the value 1. This will hold the
		running posterior probability of each possible parental genotype pair. The entries will correspond to the
		parental genotype pairs '2/2, 2/1, 2/0, 1/2, 1/1, 1/0, 0/2, 0/1, 0/0'.
	- Iterate through the sibling genotypes. For each 0, 1, or 2, multiply the nine parental posteriors by the
		posteriors calculated for 0, 1, or 2, above.
	- Renormalize the 9 probabilities by summing them together and dividing each one by the sum.
	- Initialize a 3-item array, [0,0,0], to hold the probabilities of 0, 1, and 2. For each nonzero probability
		parental genotype pair, multiply its normalized posterior probability by the likelihood of generating a
		0, 1, or 2 from that pair (which I will hardcode) as above. The final sum will give you the probability of a
		0, 1, or 2, digital sib for the observed siblings.
	"""
	geno_probs = [0, 0, 0]
	# running posterior probability of the parents being 2/2, 2/1, 2/0, 1/2, 1/1, 1/0, 0/2, 0/1, or 0/0
	total_parent_posteriors = parents_given_obs(obs_genos)
	# now compute the geno probs by multiplying the parental posteriors by the probability of generating a 0/1/2 child
	#   from that pair of parental genotypes
	for i in range(len(total_parent_posteriors)):
		geno_probs += total_parent_posteriors[i] * probs_from_parents[i]
	# now reweight geno_probs by parent posteriors
	total_parent_posteriors[total_parent_posteriors == 0] = 10**-10  # prevent numerical issues
	geno_probs[0] = np.dot(parents_given_obs(obs_genos + [0]), 1 / total_parent_posteriors)
	geno_probs[1] = np.dot(parents_given_obs(obs_genos + [1]), 1 / total_parent_posteriors)
	geno_probs[2] = np.dot(parents_given_obs(obs_genos + [2]), 1 / total_parent_posteriors)
	geno_probs = geno_probs / sum(geno_probs)
	return geno_probs


def compute_sib_geno_prob_dict(sib_genos):
	"""
	Args
	-----
	sib_genos: 3-dimensional numpy array containing genotypes of all siblings: [sib1_matrix, sib2_matrix, ... ]

	Returns
	-----
	twin_prob_dict: dict mapping observed sib genos to probability of digital twins being 0, 1, or 2
	"""
	twin_prob_dict = {}
	for i in range(len(sib_genos[0])):
		for j in range(len(sib_genos[0][i])):
			this_snp_sib_genos = sib_genos[:, i, j]  # all genotypes for this SNP for this set of sibs
			# this_snp_sib_genos = np.sort(this_snp_sib_genos)
			if str(this_snp_sib_genos) in twin_prob_dict:
				continue  # already computed probs of 0/1/2 digital sib genotype for this combo of observed genotypes
			else:
				twin_prob_dict[str(this_snp_sib_genos)] = compute_sib_geno_prob(this_snp_sib_genos)  # see above
	return twin_prob_dict


def gen_count_dict_sibs(sib_genos):
	# Helper method to count number of times each genotype combination appears in sib_genos
	# This speeds up sib generation because we can do all np random choices at once instead of one-by-one
	count_dict = {}
	num_sibs = len(sib_genos)
	for i in range(len(sib_genos[0])):
		for j in range(len(sib_genos[0][i])):
			this_snp_sib_genos = str(sib_genos[:, i, j])  # all genotypes for this SNP for this set of sibs
			if this_snp_sib_genos in count_dict:
				count_dict[this_snp_sib_genos] += num_sibs
			else:
				count_dict[this_snp_sib_genos] = num_sibs
	return count_dict


def generate_digital_sibs(sib_genos, twin_prob_dict=None, count_dict=None, seed=None, norm=False):
	"""
	Args
	-----
	sib_genos: 3-dimensional numpy array containing genotypes of all siblings: [sib1_matrix, sib2_matrix, ... ]
	twin_prob_dict: dict mapping observed sib genos to probability of digital twins being 0, 1, or 2
	count_mat: count of each sibling genotype combination, e.g. says how many (2,2) sib genotypes we have
	seed: optionally set numpy's random seed
	norm: use to normalize the returned genotype matrix of the simulated offspring
	Returns
	-----
	twin_genos numpy array with dimension equal to sib_genos, containing simulated digital sibling sets
	"""
	if seed:
		np.random.seed(seed)
	if twin_prob_dict is None:  # maps observed twin genotypes to probabilities for 0/1/2 digital sib genotypes
		twin_prob_dict = compute_sib_geno_prob_dict(sib_genos)
		count_dict = gen_count_dict_sibs(sib_genos)
	digital_sibs = np.zeros(shape=sib_genos.shape)
	num_sibs = len(sib_genos)
	geno_vals = [2, 1, 0]

	rand_draws, so_far = {}, {}
	for geno in count_dict:
		count = count_dict[geno]
		geno_probs = twin_prob_dict[geno]
		rand_draws[geno] = np.random.choice(geno_vals, p=geno_probs, size=count)
		so_far[geno] = 0

	for i in range(len(sib_genos[0])):
		for j in range(len(sib_genos[0][i])):
			this_snp_sib_genos = str(sib_genos[:, i, j])  # all genotypes for this SNP for this set of sibs
			# this_snp_sib_genos = np.sort(this_snp_sib_genos)
			# geno_probs = twin_prob_dict[str(this_snp_sib_genos)]  # probabilies of 0/1/2 genotypes for digital sibs
			for k in range(num_sibs):
				# digital_sibs[k][i][j] = np.random.choice(geno_vals, p=geno_probs)  # draw digital sib genotypes
				index = so_far[this_snp_sib_genos]
				digital_sibs[k][i][j] = rand_draws[this_snp_sib_genos][index]
				so_far[this_snp_sib_genos] += 1
	return digital_sibs, twin_prob_dict, count_dict


def shuf_haplos(genos):
	haplos = []
	for g in genos:
		if g == 0:
			haplos.append(0)
			haplos.append(0)
		elif g == 1:
			haplos.append(1)
			haplos.append(0)
		else:
			haplos.append(1)
			haplos.append(1)
	np.random.shuffle(haplos)
	new_genos = []
	for i in range(int(len(haplos) / 2)):
		new_genos.append(haplos[i*2] + haplos[i*2+1])
	return new_genos


def generate_digital_sibs_shuf(sib_genos, seed=None):
	if seed:
		np.random.seed(seed)
	digital_sibs = np.zeros(shape=sib_genos.shape)
	num_sibs = len(sib_genos)
	for i in range(len(sib_genos[0])):
		for j in range(len(sib_genos[0][i])):
			this_snp_sib_genos = sib_genos[:, i, j]  # all genotypes for this SNP for this set of sibs
			new_genos = shuf_haplos(this_snp_sib_genos)
			for k in range(num_sibs):
				digital_sibs[k][i][j] = new_genos[k]
	return digital_sibs
	

if __name__=='__main__':
	# father_dip = np.array([1,2,0,1,2,1]).reshape(1,-1)
	# mother_dip = np.array([0,1,2,2,2,1]).reshape(1,-1)
	father_dip = np.random.binomial(2, 0.3, (50000, 100))
	mother_dip = np.random.binomial(2, 0.5, (50000, 100))
	num_offspring = 1
	childgeno,parentgeno = generate_offspring(father_dip,mother_dip, num_offspring,return_parents=True, norm=True)

	print(f'childgeno shape is {childgeno.shape}')
	# print(childgeno[0])
	print(f'parents shape is {parentgeno[0].shape}')
	# print(parentgeno[0][0])
	# print(parentgeno[1][0])
