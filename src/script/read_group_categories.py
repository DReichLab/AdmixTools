import argparse
import subprocess
import sys

def is_read_group_line(line):
	return line.startswith('@RG')

def samtools_read_group_header(bam_file):
	result = subprocess.run(['samtools', 'view', '-H', bam_file], stdout=subprocess.PIPE, check=True, universal_newlines=True)
	read_group_lines = filter(is_read_group_line, result.stdout.split('\n'))
	return read_group_lines

# we are only interested in the read group id and library parts of the read group header
def parse_id_and_library(read_group_line):
	read_group_id = None
	library = None

	fields = read_group_line.split()
	for field in fields:
		if field.startswith('ID:'):
			if read_group_id:
				raise ValueError(f'multiple read group id fields {read_group_id}, {field}')
			read_group_id = field[3:]
		elif field.startswith('LB:'):
			if library:
				raise ValueError(f'multiple library id fields {library}, {field}')
			library = field[3:]
	return read_group_id, library


def samtools_read_group_ids_and_libraries(bam_file):
	return [parse_id_and_library(line) for line in samtools_read_group_header(bam_file)]

def output_lines(output_filename, lines):
	if output_filename is not None and len(output_filename) > 0:
		with open(output_filename, 'w') as f:
			for line in lines:
				print(line, file=f)
	else:
		for line in lines:
			print(line)

def read_group_categories(bam, output_filename, restrict_to=None):
	lines = []
	for read_group_id, library in samtools_read_group_ids_and_libraries(bam):
		if restrict_to is None or read_group_id in restrict_to:
			lines.append(f'{read_group_id} {len(lines)}')
	output_lines(output_filename, lines)
	return len(lines)

def read_group_categories_by_library(bam, output_filename, restrict_to=None):
	libraries = {}
	lines = []
	for read_group_id, library in samtools_read_group_ids_and_libraries(bam):
		if restrict_to is None or read_group_id in restrict_to:
			if library not in libraries:
				libraries[library] = len(libraries)
			lines.append(f'{read_group_id} {libraries[library]}')
	output_lines(output_filename, lines)
	return len(libraries)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="""Divide the readgroups for a bam into separate categories for treatment by histmake. 
	Default uses one category per read group. 
	This requires samtools.
""")

	parser.add_argument('bam', help='bam to examine read groups')
	parser.add_argument('-o', '--output_file', help='Write output to this file. If None, write to stdout')
	parser.add_argument('-l', '--library', action='store_true', help='Combine read groups for the same LB library into the same stack')
	parser.add_argument('-r', '--restrict', nargs='+', help='Use only these read groups ids')

	args = parser.parse_args()

	if args.library:
		read_group_categories_by_library(args.bam, args.output_file, args.restrict)
	else:
		read_group_categories(args.bam, args.output_file, args.restrict)
