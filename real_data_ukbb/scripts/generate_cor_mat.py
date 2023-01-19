# given a list of SNPs and a bfile, generate r^2 matrix for those SNPs with Plink
import argparse
import glob
import subprocess


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Generate r^2 matrix for given SNPs & bfile.')
    parser.add_argument('--snps', required=True, help='List of SNPs to use. Required.')
    parser.add_argument('--bfile', required=True, help='Plink bfile to calculate LD with.')
    parser.add_argument('--out', default='pruned', help='Output base name. Default: "pruned"')
    user_args = parser.parse_args()
    return user_args


if __name__ == '__main__':
    args = parseargs()
    # bfile='/u/project/sgss/UKBB/data/cal/filter4'
    subprocess.Popen(['plink', '--bfile', args.bfile, '--extract', args.snps,
                        '--make-bed', '--memory', '4096', '--out', args.out + '.TEMP_PLINK']).wait()
    subprocess.Popen(['plink', '--bfile', args.out + '.TEMP_PLINK', '--r', '--out', args.out, '--matrix']).wait()
    rm_cmd = ['rm'] + glob.glob(args.out + '.TEMP_PLINK*') + [args.out + '.log']
    subprocess.Popen(rm_cmd).wait()

