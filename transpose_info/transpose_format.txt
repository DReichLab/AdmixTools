This software package operates on file triplets:
1. genotype (geno)
2. individual (ind) This file comprises three columns: ids, sex (M/F/U), and a group label.
3. SNP (snp)

This release includes a new "transpose packed" genotype file format. The transpose packed format is very similar to the packed ancestry map format. After a short header, it contains a matrix of bit-packed individual and SNP genotype values. Transposing the matrix relative to the packed ancestry map format allows merging genotype files for multiple individuals or extracting individuals from a genotype file faster and with lower memory requirements. Additionally, the storage size of genotype file for one individual in transpose packed format is much smaller than a packed ancestry map file due to the minimum row size of 48 bytes in a packed ancestry map file.

A genotype file in packed ancestry map or transpose packed format begins with a header containing at least 48 bytes. For the transpose packed format, the header is always 48 bytes. For the packed ancestry map format, the header will be the larger of 48 bytes and the number of bytes required to store a row of SNP values.

The header contains a string with five white-space-delimited fields:
1. file type
	"GENO" packed ancestry map
	"TGENO" transpose packed
	All other values are reserved for future use.
2. number of individuals
3. number of snps
4. individual hash - This is a custom hash of the ids in the individual file.
5. snp hash - This is a custom hash of SNP ids contained in the SNP file.

The body of the genotype file is a byte-aligned matrix of SNP values for each individual. In the packed ancestry map format, each row corresponds to one SNP. In the transpose packed format, each row corresponds to one individual. Each SNP genotype occupies 2 bits for each individual. The row size of a packed ancestry map is a minimum of 48 bytes. If there is a only a single individual, only two bits of the 48 bytes, or approximately 0.5% of the file size is used.

00 homozygous reference
01 heterozygous
10 homozygous alternate
11 no data

Due to the differences in header sizes and bit padding required to byte align the rows of their respective formats, the packed ancestry map and transpose packed file sizes are slighty different for the same underlying matrix.

Building
-
The backbone of this package comes as source files. To run, you need to compile the C source files into executables.

	cd src
	make

Due to differences in the way that different Linux distributions and Mac package libraries, you may need to modify the Makefile.
Library names are specified by -l options in CFLAGS.
If libraries are stored in non-system locations, you need to specify where using the -L option.

Python script help
-
Python scripts used is this package use python3. They use argparse, and you access help by running scripts with the --help option.

python3 script/geno_header.py --help

Hashes
-
When a file triplet (geno, ind, snp) is generated, the genotype file contains hashes of the IDs in the corresponding individual file and snp file. These hashes tie the file triplet together to avoid assigning genotype values erroneously to either the wrong individual or SNP.

You can run hash checks manually using script/hash_checks.py

If you want to rename the IDs in the ind or snp file, you will need to update the corresponding hash in the genotype file. You can compute the hash for a file using
nick/nick_hash_file.py

	python3 script/nick_hash_file.py example.ind
		example.ind     4de9c1

	python3 script/nick_hash_file.py example.snp
		example.snp     f888b349

1. geno_header.py This tools allows you to inspect the header of a packed ancestry map or transpose packed genotype file and optionally to rewrite hash fields.

2. You can use convertf with option
hashcheck: NO

Transposing
-
You can convert a file between the packed ancestry map and transpose packed formats using the transpose utility or convertf.

The tranpsose utility operates purely on filenames and changes formats independent of the extension you provide. By default, this loads the entire genotype matrix into memory.
transpose input_file.geno output_file.geno

There is a low memory option that will perform the same transpose operation with a small ~1 MB memory usage but slower.
transpose --low-mem input_file.geno output_file.geno


Extracting individual(s) from the genotype file
-
To generate a subset of a large genotype file into a smaller genotype file that is easier to process:
python3 script/extract_individuals.py -i example -o out ind_1 ind_2 ind_3

This script assumes you have a triplet of files:
example.geno
example.ind
example.snp
where "ind_1", "ind_2", and "ind_3" are present in the first column of a line in example.ind.

This will create the output triplet:
out.geno
out.ind
out.snp (symlink of input snp file)

You can also use a file to pass the list of individuals:
python3 script/extract_individuals.py -i example -o out -f selected_individuals
where selected_individuals is a file containing the three lines:
ind_1
ind_2
ind_3

You should run this only on genotype files in the transpose packed format. You can convert a packed ancestry map genotype file to the transpose packed format using the above transpose utility.
