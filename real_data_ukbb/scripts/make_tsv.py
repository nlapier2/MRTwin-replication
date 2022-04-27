import argparse
import glob


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description='Make TSV containing all real data results.')
	parser.add_argument('--dir_name', required=True, help='Directory with results files. Required.')
	parser.add_argument('--outname', required=True, help='Name of output TSV file.')
	user_args = parser.parse_args()
	return user_args


def read_file(fname, exp, outc):
	res = [exp + ' --> ' + outc]
	with(open(fname, 'r')) as infile:
		line_num = 0
		mrtwin_res = ''
		for line in infile:
			line_num += 1
			if line_num < 4:  # ivw/egger/median results
				splits = line.strip().split(',')
				res.append(splits[3].split("'")[-2])
				if line_num == 1:
					mrtwin_res = splits[4].split("'")[-2]
			elif line_num == 4:
				continue
			else:  # brumpton results
				res.append(line.strip().split()[-1])
		res.append(mrtwin_res)
	return res


if __name__ == '__main__':
	args = parseargs()
	if not args.dir_name.endswith('/'):
		args.dir_name += '/'

	results = []
	all_fnames = sorted(glob.glob(args.dir_name + '*'), key=str.casefold)
	for this_fname in all_fnames:
		exposure, outcome = this_fname.split('---')
		exposure = exposure.split('/')[-1]
		outcome = outcome.split('_results.txt')[0]
		if exposure == outcome:
			continue
		results.append(read_file(this_fname, exposure, outcome))

	with(open(args.outname, 'w')) as outfile:
		outfile.write('Exposure --> Outcome\tIVW\tEgger\tMedian\tBrumpton\tMR-Twin\n')
		for this_res in results:
			outfile.write('\t'.join(this_res) + '\n')

