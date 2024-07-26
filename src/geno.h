#ifndef GENO_PACKING_H
#define GENO_PACKING_H

#include <stdbool.h>

#define GENO_HEADER_SIZE 48

// This is the position of a 2 bit unit within the byte array
typedef struct{
	size_t byte_index; // offset to find byte
	int bit_shift; // bits needed to shift [0, 2, 4, 6]
} bit_array_position;


typedef struct genotype_matrix_t genotype_matrix; // forward declaration to enable function pointer
/*
 * This a simple 2-dimensional data structure for storing genotypes by sample and SNP.
 * Each genotype is 2 bits.
 * 00
 * 01 heterozygous
 * 10
 * 11 no data
 *
 * genotypes can be stored in one of two formats
 * 1. One row = one SNP, all samples. This is the packed ancestry format, with variable length header depending on the row size (number of samples), with minumum row size 48 bytes. Each row is byte-aligned. Due to this design, a text eigenstrat file will be smaller for <48 samples.
 * 2. One row = one sample, all SNPs. This is the transpose packed format, with fixed length header of 48 bytes. This is the newer format, and is much faster for merging because the byte-aligned sample rows can simply be concatenated.
 *
 * The genotype matrix may be
 * 1. memory mapped reading from file
 * 2. allocated starting in memory
 */
typedef struct genotype_matrix_t{
	// normal has SNP rows and sample columns
	// transpose has sample rows and SNP columns
	bool transpose;
	size_t num_rows;
	size_t num_columns;
	size_t row_size_bytes; // to keep rows byte-aligned

	void (*get_position)(const genotype_matrix*, size_t sample, size_t snp, bit_array_position *);

	int individual_file_hash;
	int snp_file_hash;

	char *memory_map_start;
	size_t memory_map_length;

	unsigned char *genotypes;
} genotype_matrix;

void genotype_matrix_initialize(genotype_matrix *this);
// When not loading from a file, this allocates the genotype array in memory
void genotype_matrix_allocate(genotype_matrix *this, size_t num_samples, size_t num_snps, bool transpose);
void genotype_matrix_destructor(genotype_matrix *this);

// write the string for this matrix's header to the destination char buffer
void genotype_matrix_header(const genotype_matrix *this, bool transpose, char *dest);

// transpose: true  - write out a transpose_packed file
// transpose: false - write out a packed ancestry map file
int genotype_matrix_write_packed_file(const genotype_matrix *this, const char *filename, bool transpose);

// helper functions to locate bit array position
void row_position(size_t offset, bit_array_position *position);
void normal_position(const genotype_matrix *this, size_t sample, size_t snp, bit_array_position *position);
void transpose_position(const genotype_matrix *this, size_t sample, size_t snp, bit_array_position *position);

// Read a genotype matrix from file, either in the packed ancestry map or transpose packed format
// Optionally, check the individual and snp file hashes against the header
// the current implementation will not write changes in the matrix back out to the file
int genotype_matrix_read_file(genotype_matrix *this, const char * const genotype_filename, int flags, const int *individual_file_hash, const int *snp_file_hash);
// Read genotpye file as above in genotype_matrix_read_file, and place contents into memory
int genotype_matrix_read_file_full(genotype_matrix *this, const char * const genotype_filename, int flags, const int *individual_file_hash, const int *snp_file_hash);

int genotype_matrix_get_snp_hash(const genotype_matrix *this);
void genotype_matrix_set_snp_hash(genotype_matrix *this, int hash);
int genotype_matrix_get_individual_hash(const genotype_matrix *this);
void genotype_matrix_set_individual_hash(genotype_matrix *this, int hash);

// size of the byte-aligned genotype data, not including header
size_t genotype_matrix_byte_size(const genotype_matrix *this);
// size of the matrix
size_t genotype_matrix_num_snps(const genotype_matrix *this);
size_t genotype_matrix_num_samples(const genotype_matrix *this);

unsigned char genotype_matrix_get_genotype_internal(const genotype_matrix *this, const bit_array_position *position);
// simple accessor for 2 bit genotype data for a sample at particular SNP
unsigned char genotype_matrix_get_genotype(const genotype_matrix *this, size_t sample, size_t snp);

void array_set_genotype(unsigned char *genotypes, size_t offset, char value);
// set a genotype for the requested sample and snp
// negative values are set to 0x3
void genotype_matrix_set_genotype(genotype_matrix *this, size_t sample, size_t snp, char value);

#endif
