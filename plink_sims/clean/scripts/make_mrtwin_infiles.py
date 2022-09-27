# Given a dataframe of features and a plink assoc results, create input files for MR-Twin
import argparse
import glob
import subprocess
import numpy as np


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description='Create input files for MR-Twin given genotypes & regression results.')
	parser.add_argument('--trio_file', required=True, help='File with trio IIDs matching child to parents. Required.')
	parser.add_argument('--path_to_pheno', required=True, help='Directory with pheno files. Required.')
	parser.add_argument('--ext_bfile', required=True, help='PLINK bfile for external data. Required.')
	parser.add_argument('--child_bfile', required=True, help='PLINK bfile for child. Required.')
	parser.add_argument('--p1_bfile', required=True, help='PLINK bfile for first parent. Required.')
	parser.add_argument('--p2_bfile', required=True, help='PLINK bfile for second parent. Required.')
	parser.add_argument('--exposure_assoc', required=True, help='Plink assoc file for exposure. Required.')
	parser.add_argument('--outcome_assoc', required=True, help='Plink assoc file for outcome. Required.')
	parser.add_argument('--exposure_name', required=True, help='Outcome phenotype name. Required.')
	parser.add_argument('--outcome_name', required=True, help='Outcome phenotype name. Required.')
	parser.add_argument('--outname', required=True, help='Base name for MR-Twin input files. Required.')
	parser.add_argument('--snps', required=True, help='SNPs selected as instruments. Required.')
	user_args = parser.parse_args()
	return user_args


def read_trio_file(trio_fname):
	# track the order the trio IIDs should be in, to ensure correct child is matched up to correct parents
	child_iid_order, p1_iid_order, p2_iid_order = {}, {}, {}
	with(open(trio_fname, 'r')) as infile:
		row = 0  # row of the file
		for line in infile:
			splits = line.strip().split()
			splits = [i.split('.')[0] for i in splits]
			child, p1, p2 = splits
			child_iid_order[child] = row
			p1_iid_order[p1] = row
			p2_iid_order[p2] = row
			row += 1
	return child_iid_order, p1_iid_order, p2_iid_order


def read_famfile(fname):
	iid_dict = {}
	with(open(fname, 'r')) as infile:
		for line in infile:
			iid = line.split('\t')[1]
			iid_dict[iid] = True
	return iid_dict


def write_child_pheno(args, iid_dict):
	# using child famfile and pheno files, write pheno file for the children
	exp_vals, outc_vals = [], []  # phenotype values of exposure and outcome
	found_children = []  # list of children whose phenotype values were found
	with open(args.path_to_pheno + args.outcome_name + '.child.pheno', 'r') as infile_o, \
			open(args.path_to_pheno + args.exposure_name + '.child.pheno', 'r') as infile_e:
		o_lines = infile_o.readlines()
		e_lines = infile_e.readlines()
		eid = []
		for j in range(len(e_lines)):
			e_line = e_lines[j].strip().split(' ')
			if 'FID' in e_line:
				continue
			eid.append(int(e_line[0]))
		arr_eid = np.array(eid)
		for i in range(len(o_lines)):
			o_line = o_lines[i].strip().split(' ')
			if 'FID' in o_line:
				continue
			inds = np.where(arr_eid == int(o_line[0]))[0].tolist()
			if len(inds) == 0:
				continue
			#print(e_lines[inds[0]+1])
			print(len(e_lines))
			e_line = e_lines[inds[0]].strip().split(' ')
			eid = o_line[0]
			outcome = o_line[2]
			exposure = e_line[2]
			if eid not in iid_dict:  # check if this is one of the children in our dataset
				continue
			if exposure == '' or exposure == 'NA' or outcome == '' or outcome == 'NA':
				continue
			found_children.append(eid)
			exp_vals.append(float(exposure))
			outc_vals.append(float(outcome))
	child1_phenos = np.array([exp_vals, outc_vals])
	np.savetxt(args.outname + '.childpheno', child1_phenos)
	return found_children


def read_selected_snps(snpfile):
	selected_snps = []  # ]{}
	with(open(snpfile, 'r')) as infile:
		for line in infile:
			snp = line.strip()
			selected_snps.append(snp)
			# selected_snps[snp] = True
	return selected_snps


def select_snps_write_betas(args, selected_snps):
	beta_e, se_e, beta_o, se_o = [], [], [], []
	# select SNPs based on pval < args.pval_thresh using exposure associations
	with(open(args.exposure_assoc, 'r')) as infile:
		infile.readline()  # skip header
		for line in infile:
			splits = line.strip().split()
			snp, beta, se = splits[2], splits[8], splits[9]
			if snp not in selected_snps:
				continue
			beta, se = float(beta), float(se)
			beta_e.append(beta)
			se_e.append(se)
	# extract betas and stderrs from outcome associations for SNPs selected above
	with(open(args.outcome_assoc, 'r')) as infile:
		infile.readline()  # skip header
		for line in infile:
			splits = line.strip().split()
			snp, beta, se = splits[2], splits[8], splits[9]
			if snp not in selected_snps:
				continue
			beta, se = float(beta), float(se)
			beta_o.append(beta)
			se_o.append(se)
	betahat_eo = np.array([beta_e, se_e, beta_o, se_o])
	np.savetxt(args.outname + '.betahatEO', betahat_eo)


def get_genotypes(splits, flip):  # combine plink 12 alleles into 0/1/2 genotypes
	genos = []
	for i in range(int(len(splits) / 2)):
		genotype = 2 - (int(splits[2*i]) + int(splits[2*i+1]) - 2)
		if i in flip:  # flip genotype since effect allele in this bfile doesn't match external dataset
			genotype = 2 - genotype
		genos.append(float(genotype))
	return genos


def write_geno_file(args, bfile, iid_order, extension, snps_kept, found_children, flip_snps, pre_reorder=None):
	"""
	Write genotype files given SNPs we have selected. Slightly different for parents vs children.
	For children, we make sure only children (rows) in the pheno file are kept, and the order matches the pheno file.
	For parents, the order has already been determined by the child run, so we just use the same rows to maintain order.
	"""
	reorder = []  # ordered rows of people to keep (same order as pheno file)
	# write snps to keep, then create plink PED file with the extract command
	with(open(args.outname + '.TEMP_KEEP', 'w')) as outfile:
		for snp in snps_kept:
			outfile.write(snp + ' ' + snp + '\n')
	subprocess.Popen(['plink', '--bfile', bfile, '--extract', args.outname + '.TEMP_KEEP',
						'--recode', '12', '--memory', '4096', '--out', args.outname + '.TEMP_PLINK']).wait()
	flip_rows = find_flips_in_map(flip_snps, args.outname + '.TEMP_PLINK.map')

	# extract genotypes from ped file
	all_genos, iid_rows = [[] for i in range(len(iid_order))], {}
	# row = 0
	with(open(args.outname + '.TEMP_PLINK.ped', 'r')) as infile:
		for line in infile:
			splits = line.strip().split()
			iid = splits[0]
			ordered_row = iid_order[iid]
			if pre_reorder is None:  # if child, keep track of rows we need to keep (b/c in pheno file)
				iid_rows[iid] = ordered_row  # row
			genos = get_genotypes(splits[6:], flip_rows)
			all_genos[ordered_row] = genos  # all_genos.append(genos)
			# row += 1
	if pre_reorder is None:  # if child, reorder rows to match pheno file order
		reorder = [iid_rows[iid] for iid in found_children]  # figure out order of iids (same order as pheno file)
		all_genos = np.array(all_genos)[reorder]
	else:  # if parent, search through child files specifies reordering
		all_genos = np.array(all_genos)[pre_reorder]
	np.savetxt(args.outname + extension, all_genos)
	# clean up files
	plink_fnames = glob.glob(args.outname + '.TEMP_PLINK.*')
	subprocess.Popen(['rm', args.outname + '.TEMP_KEEP'] + plink_fnames).wait()
	return reorder


def check_flips(args, ext_bfile, other_bfile):
	flip_snps = {}
	ext_snp2a1, other_snp2a1 = {}, {}
	freqext, freqother = args.outname + '.TEMPFREQEXT', args.outname + '.TEMPFREQCHILD'
	subprocess.Popen(['plink', '--bfile', ext_bfile, '--freq', '--memory', '4096', '--out', freqext]).wait()
	subprocess.Popen(['plink', '--bfile', other_bfile, '--freq', '--memory', '4096', '--out', freqother]).wait()

	with(open(freqext + '.frq')) as infile:
		infile.readline()
		for line in infile:
			splits = line.strip().split()
			snp, a1 = splits[1], splits[2]
			ext_snp2a1[snp] = a1
	with(open(freqother + '.frq')) as infile:
		infile.readline()
		for line in infile:
			splits = line.strip().split()
			snp, a1 = splits[1], splits[2]
			other_snp2a1[snp] = a1

	for snp in other_snp2a1:
		if snp in ext_snp2a1 and other_snp2a1[snp] != ext_snp2a1[snp]:
			flip_snps[snp] = True

	# clean up files
	plink_fnames = glob.glob(freqext + '*') + glob.glob(freqother + '*')
	subprocess.Popen(['rm'] + plink_fnames).wait()
	return flip_snps


def find_flips_in_map(flip_snps, mapfile):
	flip_rows = {}
	with(open(mapfile, 'r')) as infile:
		row = 0
		for line in infile:
			snp = line.split('\t')[1]
			if snp in flip_snps:
				flip_rows[row] = True
			row += 1
	print('Flipping ' + str(len(flip_rows)) + ' SNPs.')
	return flip_rows


if __name__ == '__main__':
	main_args = parseargs()
	main_child_iid_order, main_p1_iid_order, main_p2_iid_order = read_trio_file(main_args.trio_file)
	iid_dict_child = read_famfile(main_args.child_bfile + '.fam')
	main_found_children = write_child_pheno(main_args, iid_dict_child)
	main_selected_snps = read_selected_snps(main_args.snps)
	select_snps_write_betas(main_args, main_selected_snps)

	flip_snps_child = check_flips(main_args, main_args.ext_bfile, main_args.child_bfile)
	flip_snps_p1 = check_flips(main_args, main_args.ext_bfile, main_args.p1_bfile)
	flip_snps_p2 = check_flips(main_args, main_args.ext_bfile, main_args.p2_bfile)

	main_reorder = write_geno_file(main_args, main_args.child_bfile, main_child_iid_order,
									'.childgeno', main_selected_snps, main_found_children, flip_snps_child)
	write_geno_file(main_args, main_args.p1_bfile, main_p1_iid_order,
					'.pa1geno', main_selected_snps, main_found_children, flip_snps_p1, main_reorder)
	write_geno_file(main_args, main_args.p2_bfile, main_p2_iid_order,
					'.pa2geno', main_selected_snps, main_found_children, flip_snps_p2, main_reorder)


