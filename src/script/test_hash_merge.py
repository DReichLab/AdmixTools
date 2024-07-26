"""
Testing suite for:
1. convertf between eigenstrat, packed ancestry map formats, and transpose packed formats
2. merges with merge_transpose
3. transpose (separate program) between transpose packed format and packed ancestry map
4. generating genotype subsets for single individuals

Geno/ind/snp triplets should exist in test directory.
Implicitly assumed that each test file within a triplet shares the same prefix (e.g. test1.geno, 
test1.snp, test1.ind), has the geno/ind/snp suffixes, and is stored in eigenstrat format.
Also assumed that every triplet shares the same snp set.
"""

import filecmp
import unittest
import tempfile
import os.path
import shutil
import subprocess
from itertools import zip_longest
from pathlib import Path

from convertf import convertf
from geno_header import geno_header
from hash_checks import hash_checks
import nick_hash_file
from merge_transpose import merge_geno_snp_ind
from extract_individuals import extract_single_sample_stem, extract_multiple_individuals

class TestHashMerge(unittest.TestCase):
	# all test cases use the same SNP file
	snp_hash = 'f888b349'
	num_snps = 17

	num_samples =  {
		'test1': 5,
		'test2': 7,
		'test3': 4
	}

	ind_hash = {
		'test1':	'4de9c1',
		'test2':	'69676cfd',
		'test3':	'607dc',
		'test12':	'c0b7c94c',
		'test13':	'4b6fd81d',
		'test123':	'f904a810'
	}

	GENO = 'GENO'
	TGENO = 'TGENO'
	EIGENSTRAT = 'eigenstrat'
	PACKEDANCESTRYMAP = 'packedancestrymap'
	TRANSPOSE_PACKED = 'transpose_packed'
	header_format_map = {
		PACKEDANCESTRYMAP: GENO,
		TRANSPOSE_PACKED: TGENO
	}

	executable_dir = Path(__file__).resolve().parent.parent
	test_file_directory = executable_dir / 'test'

	def genotype_file_check(self, geno_stem, expected_file_format, expected_num_samples, expected_num_snps, expected_ind_hash, expected_snp_hash):
		file_format, num_samples, num_snps, ind_hash, snp_hash = geno_header(geno_stem + '.geno')
		self.assertEqual(expected_file_format, file_format)
		self.assertEqual(expected_num_samples, num_samples)
		self.assertEqual(expected_num_snps, num_snps)
		self.assertEqual(expected_ind_hash, ind_hash)
		self.assertEqual(expected_snp_hash, snp_hash)
		if expected_file_format == self.GENO or expected_file_format == self.TGENO:
			hash_checks(geno_stem + '.geno', geno_stem + '.ind', geno_stem + '.snp')

	def _ind_file_comparison(self, x, y):
		with open(x) as file1, open(y) as file2:
			for (line1, line2) in zip_longest(file1.readlines(), file2.readlines()):
				if line1 and line2:
					fields1 = line1.strip().split()
					fields2 = line2.strip().split()
					for field1, field2 in zip_longest(fields1, fields2):
						if field1 != field2:
							return False
				else:
					return False
		return True

	def _ind_comparison(self, stem1, stem2):
		stem1_path = str(self.test_file_directory / stem1)
		ind_file1 = stem1_path + '.ind'
		stem2_path = str(self.test_file_directory / stem2)
		ind_file2 = stem2_path + '.ind'
		return self._ind_file_comparison(ind_file1, ind_file2)

	def test_ind_1_reflexive(self):
		self.assertTrue(self._ind_comparison('test1', 'test1'))
	def test_ind_2_reflexive(self):
		self.assertTrue(self._ind_comparison('test2', 'test2'))
	def test_ind_3_reflexive(self):
		self.assertTrue(self._ind_comparison('test3', 'test3'))
	def test_ind_12_reflexive(self):
		self.assertTrue(self._ind_comparison('test12', 'test12'))
	def test_ind_13_reflexive(self):
		self.assertTrue(self._ind_comparison('test13', 'test13'))
	def test_ind_123_reflexive(self):
		self.assertTrue(self._ind_comparison('test123', 'test123'))

	def test_ind_1_2(self):
		self.assertFalse(self._ind_comparison('test1', 'test2'))
	def test_ind_1_123(self):
		self.assertFalse(self._ind_comparison('test1', 'test123'))
	def test_ind_2_123(self):
		self.assertFalse(self._ind_comparison('test2', 'test123'))

	def test_eigenstrat_geno_header_fail(self):
		eigenstrat_geno_file = str(self.test_file_directory / 'test1.geno')
		self.assertRaises(ValueError, geno_header, eigenstrat_geno_file)

	# Now that convertf supports the transpose_packed format use two-way (convert to one format and back) transpose tests
	def convertf_checks(self, test_key, geno_formats):
		# end with eigenstrat format to check against inputs
		if geno_formats[-1] != self.EIGENSTRAT:
			geno_formats.append(self.EIGENSTRAT)
		with tempfile.TemporaryDirectory() as temp_directory:
			input_stem = str(self.test_file_directory / test_key)
			last_stem = input_stem
			for index, next_format in enumerate(geno_formats):
				next_stem = str(Path(temp_directory) / str(index))
				convertf(self.executable_dir, last_stem, next_stem, next_format)
				# check header contents if there is one
				if next_format in self.header_format_map:
					expected_header = self.header_format_map[next_format]
					self.genotype_file_check(next_stem, expected_header, self.num_samples[test_key], self.num_snps, self.ind_hash[test_key], self.snp_hash)
				# eigenstrat is tested against inputs
				elif next_format == self.EIGENSTRAT:
					for extension in ['.geno', '.snp']:
						self.assertTrue(filecmp.cmp(input_stem + extension, next_stem + extension, shallow=False))
						print('checking ' + extension)

	def _transpose(self, stem1, stem2, low_mem=False):
		transpose_command_array = [self.executable_dir / 'transpose']
		if low_mem:
			transpose_command_array += ['--low-mem']
		transpose_command_array += [stem1 + '.geno', stem2 + '.geno']
		subprocess.run(transpose_command_array, check=True)

	def test_hash_test1_packedancestrymap(self):
		self.convertf_checks('test1', [self.PACKEDANCESTRYMAP])

	def test_hash_test1_transpose(self):
		self.convertf_checks('test1', [self.TRANSPOSE_PACKED])

	def test_hash_test1_transpose_packedancestrymap(self):
		self.convertf_checks('test1', [self.TRANSPOSE_PACKED, self.PACKEDANCESTRYMAP])

	def test_hash_test1_packedancestrymap_transpose(self):
		self.convertf_checks('test1', [self.PACKEDANCESTRYMAP, self.TRANSPOSE_PACKED])

	def test_hash_test2_packedancestrymap(self):
		self.convertf_checks('test2', [self.PACKEDANCESTRYMAP])

	def test_hash_test2_transpose(self):
		self.convertf_checks('test2', [self.TRANSPOSE_PACKED])

	def test_hash_test2_transpose_packedancestrymap(self):
		self.convertf_checks('test2', [self.TRANSPOSE_PACKED, self.PACKEDANCESTRYMAP])

	def test_hash_test2_packedancestrymap_transpose(self):
		self.convertf_checks('test2', [self.PACKEDANCESTRYMAP, self.TRANSPOSE_PACKED])

	def test_hash_test3_packedancestrymap(self):
		self.convertf_checks('test3', [self.PACKEDANCESTRYMAP])

	def test_hash_test3_transpose(self):
		self.convertf_checks('test3', [self.TRANSPOSE_PACKED])

	def test_hash_test3_transpose_packedancestrymap(self):
		self.convertf_checks('test3', [self.TRANSPOSE_PACKED, self.PACKEDANCESTRYMAP])

	def test_hash_test3_packedancestrymap_transpose(self):
		self.convertf_checks('test3', [self.PACKEDANCESTRYMAP, self.TRANSPOSE_PACKED])

	def _convert_to_transpose_packed(self, input_stem_names, temp_directory):
		transposed_stems = []
		num_samples = 0
		for input_stem_name in input_stem_names:
			input_stem_path = str(self.test_file_directory / input_stem_name)
			transpose_stem = str(Path(temp_directory) / input_stem_name)
			convertf(self.executable_dir, input_stem_path, transpose_stem, self.TRANSPOSE_PACKED)
			transposed_stems.append(transpose_stem)
			num_samples += self.num_samples[input_stem_name]
		return transposed_stems, num_samples

	def _merge_transpose_base(self, input_stem_names, expected_stem_name, low_mem):
		expected_stem = str(self.test_file_directory / expected_stem_name)

		with tempfile.TemporaryDirectory() as temp_directory:
			# convert input eigenstrat to transpose_packed format
			transposed_stems, num_samples = self._convert_to_transpose_packed(input_stem_names, temp_directory)

			# merge
			merged_transpose_stem = str(Path(temp_directory) / 'merged_transpose')
			merge_geno_snp_ind(self.executable_dir / 'merge_transpose', transposed_stems, merged_transpose_stem)
			self.genotype_file_check(merged_transpose_stem, self.TGENO, num_samples, self.num_snps, self.ind_hash[expected_stem_name], self.snp_hash)

			# convert back to eigenstrat, first to packedancestrymap
			merged_packedancestrymap_stem = str(Path(temp_directory) / 'merged_packedancestrymap')
			self._transpose(merged_transpose_stem, merged_packedancestrymap_stem, low_mem)
			for extension in ['.ind', '.snp']:
				shutil.copy(merged_transpose_stem + extension, merged_packedancestrymap_stem + extension)

			merged_eigenstrat_stem = str(Path(temp_directory) / 'merged_eigenstrat')
			convertf(self.executable_dir, merged_packedancestrymap_stem, merged_eigenstrat_stem, self.EIGENSTRAT)

			# verify file contents
			extension = '.geno'
			self.assertTrue(filecmp.cmp(expected_stem + extension, merged_eigenstrat_stem + extension, shallow=False))

	def test_merge_transpose_12(self):
		self._merge_transpose_base(['test1', 'test2'], 'test12', False)

	def test_merge_transpose_12_low_mem(self):
		self._merge_transpose_base(['test1', 'test2'], 'test12', True)

	def test_merge_transpose_13(self):
		self._merge_transpose_base(['test1', 'test3'], 'test13', False)

	def test_merge_transpose_13_low_mem(self):
		self._merge_transpose_base(['test1', 'test3'], 'test13', True)

	def test_merge_transpose_123(self):
		self._merge_transpose_base(['test1', 'test2', 'test3'], 'test123', False)

	def test_merge_transpose_123_low_mem(self):
		self._merge_transpose_base(['test1', 'test2', 'test3'], 'test123', True)

	def _hash_ind_failures(self, input_stem_names):
		with tempfile.TemporaryDirectory() as temp_directory:
			transposed_stems, num_samples = self._convert_to_transpose_packed(input_stem_names, temp_directory)
			geno_header(transposed_stems[0] + '.geno', individual_hash='bad123')

			merged_transpose_stem = str(Path(temp_directory) / 'merged_transpose')
			self.assertRaises(ValueError, merge_geno_snp_ind, self.executable_dir / 'merge_transpose', transposed_stems, merged_transpose_stem)

	def test_merge_transpose_12_hash_fail(self):
		self._hash_ind_failures(['test1', 'test2'])

	def test_merge_transpose_13_hash_fail(self):
		self._hash_ind_failures(['test1', 'test3'])

	def test_merge_transpose_123_hash_fail(self):
		self._hash_ind_failures(['test1', 'test2', 'test3'])

	def _transpose_base(self, input_stem_name):
		with tempfile.TemporaryDirectory() as temp_directory:
			# convert input eigenstrat to transpose_packed format
			transposed_stems, num_samples = self._convert_to_transpose_packed([input_stem_name], temp_directory)
			transpose_stem = transposed_stems[0]
			packedancestrymap_stem = str(Path(temp_directory) / 'merged_packedancestrymap')
			self._transpose(transpose_stem, packedancestrymap_stem)
			transpose_stem2 = str(Path(temp_directory) / 'transpose_x2')
			self._transpose(packedancestrymap_stem, transpose_stem2)
			self.assertTrue(filecmp.cmp(transpose_stem + '.geno', transpose_stem2 + '.geno'))

	def test_transpose_x2_test1(self):
		self._transpose_base('test1')

	def test_transpose_x2_test2(self):
		self._transpose_base('test2')

	def test_transpose_x2_test3(self):
		self._transpose_base('test3')

	def _single_sample_base(self, input_stem_name, expected_stem_name, ind_identifier):
		expected_stem = str(self.test_file_directory / expected_stem_name)
		with tempfile.TemporaryDirectory() as temp_directory:
			# convert input eigenstrat to transpose_packed format
			transposed_stems, num_samples = self._convert_to_transpose_packed([input_stem_name], temp_directory)
			transposed_stem = transposed_stems[0]

			single_stem = Path(temp_directory) / 'single'
			geno_single_exec = Path(self.executable_dir) / 'geno_single'
			extract_single_sample_stem(transposed_stem, single_stem, ind_identifier, geno_single_exec)
			# transpose to packed ancestry map
			packedancestrymap_stem = str(Path(temp_directory) / 'packedancestrymap_stem')
			self._transpose(str(single_stem), packedancestrymap_stem)
			# convertf back to eigenstrat
			shutil.copy(str(self.test_file_directory / (input_stem_name + '.snp')), packedancestrymap_stem + '.snp')
			shutil.copy(str(single_stem) + '.ind', packedancestrymap_stem + '.ind')
			final_stem = str(Path(temp_directory) / 'final')
			convertf(self.executable_dir, packedancestrymap_stem, final_stem, self.EIGENSTRAT)

			self.assertTrue(self._ind_comparison(final_stem, expected_stem))
			# verify file contents
			extension = '.geno'
			self.assertTrue(filecmp.cmp(expected_stem + extension, final_stem + extension, shallow=False))

	def test_single_sample_2_H(self):
		self._single_sample_base('test2', 'H', 'H')

	def test_single_sample_2_J(self):
		self._single_sample_base('test2', 'J', 'J')

	def test_single_sample_2_K(self):
		self._single_sample_base('test2', 'K', 'K')

	def _extract_base(self, input_stem_name, individuals_list, expected_stem_name):
		with tempfile.TemporaryDirectory() as temp_directory:
			# convert to transpose packed
			input_stem = str(self.test_file_directory / input_stem_name)
			expected_stem = str(self.test_file_directory / expected_stem_name)
			transpose_stem = temp_directory + '/transpose'
			extracted_stem = temp_directory + '/extracted'
			final_stem = temp_directory + '/final'
			convertf(self.executable_dir, input_stem, transpose_stem, self.TRANSPOSE_PACKED)
			# extract samples
			geno_single = self.executable_dir / 'geno_single'
			merge_transpose = self.executable_dir / 'merge_transpose'
			extract_multiple_individuals(transpose_stem, extracted_stem, individuals_list, geno_single, merge_transpose)
			# convert back to eigenstrat
			convertf(self.executable_dir, extracted_stem, final_stem, self.EIGENSTRAT)
			# check
			self._ind_comparison(final_stem, expected_stem)
			extension= '.geno'
			self.assertTrue(filecmp.cmp(final_stem + extension, expected_stem + extension, shallow=False))

	def test_extract_1(self):
		self._extract_base('test123', ['A', 'B', 'C', 'D', 'E'], 'test1')

	def test_extract_2(self):
		self._extract_base('test123', ['F', 'G', 'H', 'I', 'J', 'K', 'L'], 'test2')

	def test_extract_3(self):
		self._extract_base('test123', ['M', 'N', 'O', 'P'], 'test3')

	# This is a test using convertf ignore feature
	def ignore_base(self, geno_format, working_directory = None):
		with tempfile.TemporaryDirectory() as temp_directory:
			input_stem = str(Path(self.test_file_directory) / 'ignore' / 'start')
			expected_stem = str(Path(self.test_file_directory) / 'ignore' / 'expected')

			directory = working_directory if working_directory else temp_directory
			ignore_removed_stem = directory + '/ignore_removed'
			finished_stem = directory + '/finished'

			convertf(self.executable_dir, input_stem, ignore_removed_stem, geno_format, par_directory=working_directory)
			expected_header = self.header_format_map[geno_format]
			self.genotype_file_check(ignore_removed_stem, expected_header, 2, 17, '482', self.snp_hash)
			convertf(self.executable_dir, ignore_removed_stem, finished_stem, self.EIGENSTRAT)
			for extension in ['.geno', '.ind']:
				self.assertTrue(filecmp.cmp(finished_stem + extension, expected_stem + extension, shallow=False))
				print('checking ' + extension)

	def test_convertf_ignore_packedancestrymap(self):
		self.ignore_base(self.PACKEDANCESTRYMAP)

	def test_convertf_ignore_transpose(self):
		self.ignore_base(self.TRANSPOSE_PACKED)

if __name__ == '__main__':
	unittest.main()
