import argparse
import subprocess
from pathlib import Path

inferred_executable = Path(__file__).resolve().parent.parent / 'nickhash'

# Nick's implementation hashes the first column and ignores everything else
def hash_single_file(filename, executable=inferred_executable):
	result = subprocess.run([executable, '-i', filename], check=True, stdout=subprocess.PIPE, universal_newlines=True)
	hash_str = result.stdout.strip()
	hash_value = int(hash_str, 16)
	return f'{hash_value:x}'

def ind_hash(filename, executable=inferred_executable):
	return hash_single_file(filename, executable)

def snp_hash(filename, executable=inferred_executable):
	return hash_single_file(filename, executable)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Compute Nick's custom hash on a ind or snp file.")
	parser.add_argument('--exe', help=f"path to nickhash executable", default=inferred_executable)
	parser.add_argument('files', help="Individual (ind) file to compute hash", nargs='+')

	args = parser.parse_args()

	for f in args.files:
		hash_value = hash_single_file(f, args.exe)
		print(f'{f}\t{hash_value}')
