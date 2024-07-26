import argparse
import sys
from geno_header import geno_header
from pathlib import Path
from nick_hash_file import ind_hash, snp_hash

# returns True if hash checks pass
def hash_checks(geno_file, ind_file, snp_file):
	individual_file_hash = ind_hash(ind_file)
	snp_file_hash = snp_hash(snp_file)

	file_format, num_samples, num_snps, geno_ind_hash, geno_snp_hash = geno_header(geno_file)
	mismatch = False
	if individual_file_hash != geno_ind_hash:
		raise ValueError(f'individual hash mismatch {individual_file_hash} {geno_ind_hash}')
	if snp_file_hash != geno_snp_hash:
		raise ValueError(f'snp hash mismatch {snp_file_hash} {geno_snp_hash}')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Inspect ind and snp hashes for genotype file.")

	parser.add_argument('geno', help='Genotype file for hash checks. File should be in packed ancestry map or transpose packed format.')
	parser.add_argument('-i', '--individual_file', help='Only needed if file name stem is different or filename does not end in .ind')
	parser.add_argument('-s', '--snp_file', help='Only needed if file name stem is different or filename does not end in .snp')

	args = parser.parse_args()

	geno_path = Path(args.geno)
	parent_path = geno_path.resolve().parent

	individual_file = args.individual_file if args.individual_file else (parent_path / (geno_path.stem + '.ind'))
	snp_file = args.snp_file if args.snp_file else (parent_path / (geno_path.stem + '.snp'))

	hash_checks(args.geno, individual_file, snp_file)
	print('hashes match')
