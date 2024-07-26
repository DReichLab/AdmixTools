import argparse
from pathlib import Path
import subprocess
import sys

def hash_program(f):
	hash_function = 'sha256sum'
	result = subprocess.run([hash_function, f], check=True, universal_newlines=True, stdout=subprocess.PIPE) # for python 3.7, universal_newlines -> text
	hash_value, filename = result.stdout.strip().split()
	return hash_value

def link_identical_files(reference, to_compare, verbose=False):
	reference_hash = hash_program(reference)
	reference_path = Path(reference)
	if verbose:
		print(f'reference hash {reference_hash}', file=sys.stderr)
	for f in to_compare:
		comparison_hash = hash_program(f)
		if verbose:
			print(f'comparison hash {comparison_hash}', file=sys.stderr)
		if reference_hash == comparison_hash:
			print(f'symlinking {f}', file=sys.stderr)
			path = Path(f)
			path.unlink()
			path.symlink_to(reference_path)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Compare file(s) to a reference, and if the same, then symbolic link compared files. This is originally designed for snp files, which are output upon every pulldown")

	parser.add_argument('-r', '--ref', required=True, help="Reference file")
	parser.add_argument('-v', '--verbose', action='store_true', help="verbose: print hashes")
	parser.add_argument('files', nargs='*', help='Files to compare to reference')

	args = parser.parse_args()

	link_identical_files(args.ref, args.files, args.verbose)
