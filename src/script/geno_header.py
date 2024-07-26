import argparse
from pathlib import Path
import re

def geno_header(filename, individual_hash=None, snp_hash=None, header_size=48, **kwargs):
	mode = 'rb'
	rewrite_header = individual_hash is not None or snp_hash is not None
	if rewrite_header:
		mode = 'rb+'
	with open(filename, mode) as f:
		original_header = f.read(header_size)

		truncated_header = original_header.decode('utf-8').strip('\x00')
		fields = truncated_header.split()
		file_format = fields[0]
		if file_format != 'GENO' and file_format != 'TGENO':
			raise ValueError(f'Unsupported file format {file_format}. Only GENO and TGENO are supported.')

		num_samples = int(fields[1])
		num_snps = int(fields[2])
		individual_hash_int = int(fields[3], 16)
		snp_hash_int = int(fields[4], 16)

		print(truncated_header)
		#print(len(truncated_header))

		if rewrite_header:
			if individual_hash is not None:
				individual_hash_int = int(individual_hash, 16)
			if snp_hash is not None:
				snp_hash_int = int(snp_hash, 16)

			header_string = f'{file_format} {num_samples:7d} {num_snps:7d} {individual_hash_int:x} {snp_hash_int:x}'

			print(header_string)
			header = bytes(header_string, 'utf-8')
			bytes_to_pad = header_size - len(header)
			if bytes_to_pad < 0:
				raise ValueError(f'Required header size is too big {len(header)} for fixed size header {header_size}')
			else:
				f.seek(0)
				f.write(header)
				# pad with null
				if bytes_to_pad > 0:
					f.write(bytes('\x00' * bytes_to_pad, 'utf-8'))
	return file_format, num_samples, num_snps, f'{individual_hash_int:x}', f'{snp_hash_int:x}'

def add_to_kwargs(kwargs, args, property_name):
	try:
		value = getattr(args, property_name)
		kwargs[property_name] = value
	except AttributeError:
		pass

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Inspect a genotype header file. Optionally replace the individual or SNP file hashes in place.")

	parser.add_argument('geno', help='Genotype file to have its header inspected or overwritten. File should be in packed ancestry map or transpose packed format.')
	parser.add_argument('-i', '--individual_hash', help="hex string to replace Nick's 32 bit hash of the individuals from the ind file")
	parser.add_argument('-s', '--snp_hash', help="hex string to replace Nick's 32 bit hash of the snp file")
	parser.add_argument('-n', '--header_size', type=int, default=48, help="size of header to overwrite")

	args = parser.parse_args()

	kwargs = {}
	for argument_name in ['individual_hash', 'snp_hash', 'header_size']:
		add_to_kwargs(kwargs, args, argument_name)

	geno_header(args.geno, **kwargs)
