import argparse
import subprocess
import tempfile

from hash_checks import hash_checks
from geno_header import geno_header
from merge_transpose import transpose_packed_merge, merge_ind_files
from nick_hash_file import ind_hash, snp_hash
from pathlib import Path

# extract single individual from ind_filename and put it into output_filename
def individual_index(ind_filename, output_filename, single_id):
	with open(ind_filename) as f:
		for index, line in enumerate(f):
			identifier, sex, group = line.strip().split()
			if single_id == identifier:
				with open(output_filename, 'w') as out:
					print(line, end='', file=out)
				return index
		raise ValueError(f'{single_id} not found in individuals file')

def extract_single_sample_stem(input_stem, output_stem, single_id, geno_single_exec, input_geno_extension='geno'):
	input_stem_str = str(input_stem)
	output_stem_str = str(output_stem)
	input_ind = input_stem_str + '.ind'
	input_geno = input_stem_str + '.' + input_geno_extension
	print(input_geno) 
	output_ind = output_stem_str + '.ind'
	output_geno = output_stem_str + '.geno'

	index = individual_index(input_ind, output_ind, single_id)
	print(index)
	individuals_hash = ind_hash(output_ind)
	subprocess.run([str(geno_single_exec), input_geno, str(index), output_geno], check=True)
	geno_header(output_geno, individual_hash=individuals_hash)

def extract_multiple_individuals(input_stem, output_stem, individuals_list, geno_single_exec, merge_transpose_exec, input_geno_extension='geno', working_directory=None):
	with tempfile.TemporaryDirectory() as temp_directory:
		directory = working_directory if working_directory is not None else temp_directory

		# extract single individuals to temporary files
		component_stems = str(Path(temp_directory) / "component_stems")
		with open(component_stems, 'w') as f:
			for index, individual_id in enumerate(individuals_list):
				component_stem = f'{directory}/{index}'
				extract_single_sample_stem(input_stem, component_stem, individual_id, geno_single_exec, input_geno_extension)
				print(component_stem, file=f)
		# merge temporary_files
		# ind
		output_individual_file = output_stem + '.ind'
		merge_ind_files(output_individual_file, component_stems)
		individual_hash_str = ind_hash(output_individual_file)
		# snp
		input_snp_file = input_stem + '.snp'
		snp_hash_str = snp_hash(input_snp_file)
		output_snp_filename = output_stem + '.snp'
		output_snp_file = Path(output_snp_filename)
		if output_snp_file.exists():
			if not output_snp_file.samefile(input_snp_file):
				raise FileExistsError(output_snp_file)
		else:
			input_snp_file_fullpath = Path(input_snp_file).resolve()
			output_snp_file.symlink_to(input_snp_file_fullpath)
		# geno
		transpose_packed_merge(merge_transpose_exec, component_stems, output_stem, ind_file_hash_str=individual_hash_str, snp_file_hash_str=snp_hash_str)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract individuals from a genotype file. This will be fast from a transpose packed formatted file. Input is triplet of (.geno, .ind, .snp) files with the same file stem. Output is a pair of .geno and .ind files with the same stem and a symlinked .snp to complete the triplet.")

	parser.add_argument('ids', help='Individual identifiers from .ind file', nargs='*')
	parser.add_argument('-f', '--file_list', help='Individual identifiers from .ind file, loaded from a file.')
	parser.add_argument('-i', '--input_stem', help='Stem of genotype and ind file', required=True)
	parser.add_argument('-o', '--output_stem', help='Stem of genotype and ind file', required=True)
	parser.add_argument('-e', '--geno_single_exec', help='geno_single executable', default=Path(__file__).resolve().parent.parent / 'geno_single') # executable is up one level
	parser.add_argument('-m', '--merge_transpose_exec', help='merge_transpose executable', default=Path(__file__).resolve().parent.parent / 'merge_transpose') # executable is up one level
	parser.add_argument('-x', '--hash_ignore', action='store_true', help='skip ind and snp hash checks')
	parser.add_argument('-t', '--tgeno', action='store_true', help='use tgeno extension instead of geno for input')
	parser.add_argument('--temp_directory', help='Use this directory for temporary output files. This is useful if you want to save files for inspection.')

	args = parser.parse_args()

	individuals_list = args.ids

	if args.file_list:
		with open(args.file_list) as f:
			for line in f:
				individual_id = line.strip().split()[0]
				if len(individual_id) > 0 and not individual_id.startswith('#'):
					individuals_list += [individual_id]

	input_geno_extension = 'tgeno' if args.tgeno else 'geno'
	if not args.hash_ignore:
		hash_checks(args.input_stem + '.' + input_geno_extension, args.input_stem + '.ind', args.input_stem + '.snp')
	extract_multiple_individuals(args.input_stem, args.output_stem, individuals_list, args.geno_single_exec, args.merge_transpose_exec, input_geno_extension, args.temp_directory)
