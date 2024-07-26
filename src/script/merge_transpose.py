import argparse
import filecmp
import fileinput
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from nick_hash_file import ind_hash, snp_hash
from pathlib import Path
from geno_header import geno_header

# return true if snp files are the same
# return false otherwise
def compare_snp_files(snp_file_array):
	for other_file in snp_file_array[1:]:
		if filecmp.cmp(snp_file_array[0], other_file, shallow=False) == False:
			return False
	return True

def check_snp_file(snp_filename):
	BASES = set('ACGT')
	count = 0
	with open(snp_filename, 'r') as snps:
		for snp_line in snps:
			fields = snp_line.split()
			if len(fields) != 6:
				raise ValueError("does not have expected 6 fields: '{}'".format(snp_line))
			chromosome = fields[1]
			if not chromosome.isdigit() or (int(chromosome) < 1 or int(chromosome) > 24):
				raise ValueError("bad chromosome '{}', {}".format(snp_line, chromosome))
			
			try:
				cm = float(fields[2])
			except:
				raise ValueError("nonfloat '{}': {}".format(snp_line, fields[2]))
			
			position = fields[3]
			if not position.isdigit():
				raise ValueError("bad position '{}': {}".format(snp_line, position))
			
			allele1 = fields[4]
			allele2 = fields[5]
			if (allele1 not in BASES) or (allele2 not in BASES) or (allele1 == allele2):
				raise ValueError("bad alleles '{}': {} {}".format(snp_line, allele1, allele2))
			count += 1
	return count

# combine the contents of multiple individual files into one, starting from a file list
def merge_ind_files(output_filename, input_stem_file, max_overlap=1):
	with open(input_stem_file) as f:
		input_stems = f.readlines()
		merge_ind_files_list(output_filename, input_stems, max_overlap)

# combine the contents of multiple individual files into one, starting from a python list
def merge_ind_files_list(output_filename, input_stems, max_overlap=1):
	ind_counts = [] # individual counts for each individual file
	individuals_geno_offsets_dict = OrderedDict() # mapping of geno sets to individuals by geno offset, treating geno data as one array
	VALID_SEX = frozenset('MFU')
	with open(output_filename, 'w') as out:
		overall_line_count = 0
		for input_stem in input_stems:
			ind_filename = input_stem.strip() + '.ind'
			prior_files_total = 0 if len(ind_counts) == 0 else overall_line_count
			with open(ind_filename, 'r') as ind_file:
				for line in ind_file:
					fields = line.split()
					if len(fields) > 0:
						individual = fields[0]
						if individual not in individuals_geno_offsets_dict:
							out.write(line)
							individuals_geno_offsets_dict[individual] = [overall_line_count]
						else:
							# make sure there are no duplicate entries for an individual
							if len(individuals_geno_offsets_dict[individual]) >= max_overlap:
								raise ValueError('individual {} appears too many times'.format(individual))
							else:
								individuals_geno_offsets_dict[individual].append(overall_line_count)
						overall_line_count += 1
						# sex is one of M,F,U
						sex = fields[1]
						if sex not in VALID_SEX:
							raise ValueError('invalid sex {} in {}'.format(sex, line))

				ind_counts.append(overall_line_count - prior_files_total)
	individuals_geno_offsets = [individuals_geno_offsets_dict[x] for x in individuals_geno_offsets_dict]
	return ind_counts, individuals_geno_offsets

def transpose_packed_merge(merge_transpose, input_stem_file, output_stem, ind_file_hash_str = None, snp_file_hash_str = None, **kwargs):
	command_args = [merge_transpose, '--input', input_stem_file, '--output', output_stem + '.geno']
	if ind_file_hash_str is not None:
		command_args.append('--ind_hash')
		command_args.append(ind_file_hash_str)
	if snp_file_hash_str is not None:
		command_args.append('--snp_hash')
		command_args.append(snp_file_hash_str)
	subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

def copy_snp_file(input_stem, output_stem):
	snp_filename = input_stem + '.snp'
	num_snps = check_snp_file(snp_filename)

	output_filename = f"{output_stem}.snp"
	shutil.copyfile(snp_filename, output_filename)
	return output_filename

# Perform snp file comparison (currently SNP files must be identical)
# Merge ind and genotype files
def merge_geno_snp_ind(merge_transpose, input_stems, output_stem):
	with tempfile.TemporaryDirectory() as temp_directory:
		input_stem_filename = str(Path(temp_directory) / "input_stems")
		with open(input_stem_filename, 'w') as f:
			for input_stem in input_stems:
				print(str(Path(input_stem.strip()).resolve()), file=f)
		merge_geno_snp_ind_file(merge_transpose, input_stem_filename, output_stem)

def merge_geno_snp_ind_file(merge_transpose, input_stem_file, output_stem):
	# snp
	with open(input_stem_file) as f:
		first_input_stem = f.readline().strip()
	snp_filename = copy_snp_file(first_input_stem, output_stem)
	snp_hash_str = snp_hash(snp_filename)

	# ind
	ind_filename = output_stem + '.ind'
	merge_ind_files(ind_filename, input_stem_file)
	ind_hash_str = ind_hash(ind_filename)

	# geno
	transpose_packed_merge(merge_transpose, input_stem_file, output_stem, ind_file_hash_str=ind_hash_str, snp_file_hash_str=snp_hash_str)

	print('Merge output files written. Continuing to check SNP files.', file=sys.stderr)
	with open(input_stem_file) as f:
		snp_filenames = [f'{line.strip()}.snp' for line in f]
	if not compare_snp_files(snp_filenames):
		raise ValueError('SNP file mismatch')
	print('SNP file checks successful.', file=sys.stderr)

	# hash checks between ind, snp files and geno header
	with open(input_stem_file) as f:
		for line in f:
			input_stem = line.strip()
			ind_hash_str = ind_hash(input_stem + '.ind')
			file_format, num_samples, num_snps, geno_ind_hash, geno_snp_hash = geno_header(input_stem + '.geno')
			if geno_ind_hash != ind_hash_str:
				raise ValueError(f'{input_stem} ind hash mismatch')
			if geno_snp_hash != snp_hash_str:
				raise ValueError(f'{input_stem} snp hash mismatch')
	print('ind, SNP file hash checks successful.', file=sys.stderr)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Merge transpose_packed results with non-overlapping individuals. SNP sets must be the same.")

	parser.add_argument('--merge_transpose_executable', help='binary executable for merge_transpose', default=(Path(__file__).resolve().parent.parent / 'merge_transpose'))
	
	input_group = parser.add_mutually_exclusive_group(required=True)
	input_group.add_argument('-i', "--input_stems", help="List of stems, with each stem having three files, for example: a0.{geno,snp,ind} or b0.b1.{geno,snp,ind}.", nargs='+')
	input_group.add_argument('-f', "--file_input", help="Read from given file the List of stems, with each stem having three files, for example: a0.{geno,snp,ind} or b0.b1.{geno,snp,ind}.")

	parser.add_argument('-o', "--output_stem", help="Stem for output files: output.{geno,snp,ind}.", default='out')

	args = parser.parse_args()

	output_stem = args.output_stem

	if args.file_input:
		merge_geno_snp_ind_file(args.merge_transpose_executable, args.file_input, args.output_stem)
	elif args.input_stems:
		merge_geno_snp_ind(args.merge_transpose_executable, args.input_stems, args.output_stem)
