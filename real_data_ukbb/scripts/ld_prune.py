import argparse
import glob
import os
import subprocess


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Prune plink regression results, prioritizing most significant SNPs.')
    parser.add_argument('--gwas', required=True, help='Plink linear regression GWAS results file. Required.')
    parser.add_argument('--bfile', required=True, help='Plink bfile to calculate LD with.')
    parser.add_argument('--gwas_thresh', type=float, default=5*10**-8, help='Threshold for significant GWAS results.')
    parser.add_argument('--prune_thresh', type=float, default=0.5, help='r^2 threshold for pruning.')
    parser.add_argument('--out', default='pruned', help='Output base name. Default: "pruned"')
    parser.add_argument('--temp_basename', default='',
                        help='A basename to prepend to tempfiles, useful for parallel runs.')
    user_args = parser.parse_args()
    return user_args


def read_gwas(fname, gwas_thresh):
    # read in significant SNPs from GWAS
    chrom2snps, full_snpslist = {}, []
    with(open(fname, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            splits = line.strip().split('\t')
            chrom, snpid, pval = splits[0], splits[2], splits[11]
            if pval == 'NA':
                continue
            pval = float(pval)
            if pval > gwas_thresh:
                continue
            if chrom not in chrom2snps:
                chrom2snps[chrom] = []
            chrom2snps[chrom].append([snpid, pval])
            full_snpslist.append(snpid)
    return chrom2snps, full_snpslist


def extract_snps(bfile, snplist, snpfile='TEMP_SNPLIST', plink_out='TEMP_BFILE'):
    with(open(snpfile, 'w')) as outfile:
        for snp in snplist:
            outfile.write(snp + '\n')
    subprocess.Popen(['plink', '--bfile', bfile, '--extract', snpfile, '--make-bed', '--out', plink_out]).wait()
    subprocess.Popen(['rm', snpfile]).wait()


def read_ld_pairs(fname):
    pairs = {}
    if not os.path.exists(fname):
        return pairs
    with(open(fname, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            splits = line.strip().split()
            snp1, snp2 = splits[2], splits[5]
            if snp1 not in pairs:
                pairs[snp1] = []
            if snp2 not in pairs:
                pairs[snp2] = []
            pairs[snp1].append(snp2)
            pairs[snp2].append(snp1)
    return pairs


def make_in_out_lists(snps, pairs, inlist, outlist):
    # determine which snps to prune (out) or keep (in). adds on to existing in_list and out_list.
    snps.sort(key=lambda x: x[0])  # sort by descending pvalue
    for snp, pval in snps:
        if snp in outlist:  # don't look at snps we've already pruned
            continue
        if snp not in pairs:  # snp has no high-ld pairs
            inlist[snp] = True
            continue
        high_ld_snps = pairs[snp]  # all snps in high ld with current snp
        for pair_snp in high_ld_snps:  # prune all snps in high ld with this snp...
            if pair_snp not in inlist:  # ...unless they've already been retained due to stronger pval than this snp
                outlist[pair_snp] = True
        inlist[snp] = True
    return inlist, outlist


def write_in_out(basename, inlist, outlist):
    with(open(basename + '.in', 'w')) as outfile:
        for snp in inlist:
            outfile.write(snp + '\n')
    with(open(basename + '.out', 'w')) as outfile:
        for snp in outlist:
            outfile.write(snp + '\n')


if __name__ == '__main__':
    args = parseargs()
    in_list, out_list = {}, {}
    snpfile = args.temp_basename + 'TEMP_SNPLIST'
    tempname = args.temp_basename + 'TEMP_BFILE'
    tempnamechr = args.temp_basename + 'TEMP_BFILE_CHR'
    r2_fname = args.temp_basename + 'TEMP_R2'

    chrom2snp, full_snplist = read_gwas(args.gwas, args.gwas_thresh)
    extract_snps(args.bfile, full_snplist, snpfile=snpfile, plink_out=tempname)

    for chrom in chrom2snp:
        chrom_snps = chrom2snp[chrom]
        snps_only = [i[0] for i in chrom_snps]
        extract_snps(tempname, snps_only, snpfile=snpfile, plink_out=tempnamechr)
        subprocess.Popen(['plink', '--bfile', tempnamechr, '--r2', '--out', r2_fname, '--ld-window', '99999999',
                          '--ld-window-kb', '99999999', '--ld-window-r2', str(args.prune_thresh)]).wait()
        high_ld_pairs = read_ld_pairs(r2_fname + '.ld')
        in_list, out_list = make_in_out_lists(chrom_snps, high_ld_pairs, in_list, out_list)
        rm_cmd = ['rm'] + glob.glob(tempnamechr + '.*') + glob.glob(r2_fname + '.*')
        subprocess.Popen(rm_cmd).wait()

    write_in_out(args.out, in_list, out_list)
    rm_cmd = ['rm'] + glob.glob(tempname + '.*')
    subprocess.Popen(rm_cmd).wait()

