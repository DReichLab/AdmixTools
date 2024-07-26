import argparse
import subprocess
import tempfile
from pathlib import Path

def buildParFile(fileName, _input, _output, _format, hashCheck = False, **options):
	"""Builds parameter file for convertf"""
	parFileText = f'''genotypename:    {_input}.geno
		              snpname:         {_input}.snp
		              indivname:       {_input}.ind
		              outputformat:    {_format}
		              genotypeoutname: {_output}.geno
		              snpoutname:      {_output}.snp
		              indivoutname:    {_output}.ind'''
	if not hashCheck:
		parFileText += f'\nhashcheck:       NO'
	for option, value in options.items():
		parFileText += f'\n{option}: {value}'
	parFileText += '\n'

	with open(fileName, 'w') as f:
		f.write(parFileText)

# par_directory gives the option of using a permanent directory instead of a temporary one
# This is if you want to save the par file.
def convertf(executable_directory, input_stem, output_stem, geno_format, par_directory=None, **options):
	executable = Path(executable_directory) / 'convertf'
	with tempfile.TemporaryDirectory() as temp_directory:
		directory = par_directory if par_directory else temp_directory
		temp_par = str(Path(directory) / 'convertf.par')
		buildParFile(temp_par, input_stem, output_stem, geno_format, True, **options)
		subprocess.run([str(executable), '-p', temp_par], check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Wrapper of Nick's convertf utility")

	parser.add_argument('--convertf_dir', help="directory containing convertf", default=Path(__file__).resolve().parent.parent) # executable is up one level
	parser.add_argument('-i', '--input_stem', help="Input stem for geno,ind,snp files", required=True)
	parser.add_argument('-o', '--output_stem', help="Output stem for geno,ind,snp files", required=True)
	parser.add_argument('-f', '--geno_format', help="Format for output genotype file", choices=['eigenstrat', 'packedancestrymap', 'transpose_packed', 'packedped'], default='transpose_packed')
	parser.add_argument('-d', '--directory', help="par file working directory", default=None)
	parser.add_argument('options', nargs='*')

	args = parser.parse_args()

	keys = args.options[::2]
	values = args.options[1::2]
	if len(keys) != len(values):
		raise ValueError('additional options need to be key-value pairs')
	option_dictionary = dict(zip(keys, values))

	convertf(args.convertf_dir, args.input_stem, args.output_stem, args.geno_format, args.directory, **option_dictionary)
