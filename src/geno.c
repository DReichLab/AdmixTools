#include <assert.h>
#include <fcntl.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "geno.h"

// each entry is 2 bits
size_t padded_bytes_needed(size_t n){
	return (size_t) ceil(n/4.0);
}

void genotype_matrix_fail( const char * message){
	fprintf(stderr, "genotype matrix failure: %s\n", message);
	exit(-1);
}

void genotype_matrix_initialize(genotype_matrix *this){
	this->num_rows = 0;
	this->num_columns = 0;
	this->row_size_bytes = 0;

	this->individual_file_hash = 0;
	this->snp_file_hash = 0;

	this->memory_map_start = NULL;
	this->memory_map_length = 0;

	this->genotypes = NULL;
	this->transpose = false;
	this->get_position = 0;
}

void genotype_matrix_setup_dimensions(genotype_matrix *this, size_t num_samples, size_t num_snps, bool transpose){
	this->transpose = transpose;
	if(this->transpose){
		this->num_rows = num_samples;
		this->num_columns = num_snps;
		this->get_position = transpose_position;
		this->row_size_bytes = padded_bytes_needed(this->num_columns);
	}
	else{
		this->num_rows = num_snps;
		this->num_columns = num_samples;
		this->get_position = normal_position;
		size_t needed_bytes = padded_bytes_needed(this->num_columns);
		this->row_size_bytes = (needed_bytes > GENO_HEADER_SIZE) ? needed_bytes : GENO_HEADER_SIZE;
	}
}

void genotype_matrix_allocate(genotype_matrix *this, size_t num_samples, size_t num_snps, bool transpose){
	genotype_matrix_setup_dimensions(this, num_samples, num_snps, transpose);

	size_t size = genotype_matrix_byte_size(this);
	this->genotypes = (unsigned char *) malloc(size);
	// missing data is all 1's
	memset(this->genotypes, 0xFF, size);
	// Nick writes unused end-of-row padding to 0. We replicate that here.
	size_t max = this->row_size_bytes * 4;
	for (size_t row = 0; row < this->num_rows; row++){
		unsigned char *row_start = this->genotypes + (row * this->row_size_bytes);
		assert(max - this->num_columns < 4);
		for (size_t offset = this->num_columns; offset < max; offset++){
			array_set_genotype(row_start, offset, 0);
		}
	}
}

void genotype_matrix_destructor(genotype_matrix *this){
	if(this->memory_map_start){
		munmap(this->memory_map_start, this->memory_map_length);
	}
	else if(this->genotypes != NULL){
		free(this->genotypes);
	}

	genotype_matrix_initialize(this);
}

size_t genotype_matrix_byte_size(const genotype_matrix *this){
	return this->num_rows * this->row_size_bytes * sizeof(unsigned char);
}

int genotype_matrix_read_file(genotype_matrix *this, const char * const genotype_filename, int flags, const int  *individual_file_hash, const int *snp_file_hash){
	genotype_matrix_initialize(this);
	struct stat sb;
	size_t num_samples, num_snps;

	int fd = open(genotype_filename, O_RDONLY);
	if(fd < 0){
		return -1;
	}

	char header[GENO_HEADER_SIZE+1];
	header[GENO_HEADER_SIZE] = 0;
	char matrix_type[GENO_HEADER_SIZE+1];

	fstat(fd, &sb);
	this->memory_map_length = sb.st_size;
	this->memory_map_start = mmap(NULL, this->memory_map_length, PROT_READ, MAP_PRIVATE | flags, fd, 0);
	if(this->memory_map_start == MAP_FAILED){
		return -2;
	}
	close(fd);
	// read header
	strncpy(header, this->memory_map_start, GENO_HEADER_SIZE);
	//printf("%s\n", header);
	//printf("%d\n", strlen(header));

	sscanf(header, "%s %zd %zd %x %x", matrix_type, &num_samples, &num_snps, &(this->individual_file_hash), &(this->snp_file_hash));
	printf("matrix type is %s\n", matrix_type);
	//printf("samples %zd\n", num_samples);
	//printf("snps %zd\n", num_snps);
	bool transpose;
	if(strncmp(matrix_type, "GENO", 4) == 0){
		transpose = false;
	} else if (strncmp(matrix_type, "TGENO", 5) == 0){
		transpose = true;
	} else{
		fprintf(stderr, "Unhandled file type");
		return -3;
	}
	genotype_matrix_setup_dimensions(this, num_samples, num_snps, transpose);
	// set genotype pointer
	size_t header_size = GENO_HEADER_SIZE;
	if(!this->transpose){
		header_size = this->row_size_bytes;
	}
	this->genotypes = (unsigned char *) this->memory_map_start + header_size;
	// verify size
	size_t expected_size = header_size + genotype_matrix_byte_size(this);
	if (this->memory_map_length != expected_size){
		fprintf(stderr, "Unexpected %s size. Expected %ld was %ld\n", matrix_type, expected_size, this->memory_map_length);
		return -4;
	}

	// verify hashes
	if(individual_file_hash && (*individual_file_hash != this->individual_file_hash)){
		fprintf(stderr, "individual file hash mismatch %x %x", *individual_file_hash, this->individual_file_hash);
		return -5;
	}
	if(snp_file_hash && (*snp_file_hash != this->snp_file_hash)){
		fprintf(stderr, "snp file hash mismatch %x %x", *snp_file_hash, this->snp_file_hash);
		return -6;
	}

	return 0;
}

int genotype_matrix_read_file_full(genotype_matrix *this, const char * const genotype_filename, int flags, const int  *individual_file_hash, const int *snp_file_hash){
	int result = genotype_matrix_read_file(this, genotype_filename, flags, individual_file_hash, snp_file_hash);
	if (result < 0){
		return result;
	}	
	size_t genotype_size = genotype_matrix_byte_size(this);
	unsigned char *in_memory = malloc(genotype_size);
	memcpy(in_memory, this->genotypes, genotype_size);
	munmap(this->memory_map_start, this->memory_map_length);
	this->memory_map_start = NULL;
	this->genotypes = in_memory;
	return 0;
}

void genotype_matrix_header(const genotype_matrix *this, bool transpose, char *dest){
	char *matrix_type = transpose ? "TGENO" : "GENO";
	size_t num_samples = this->transpose ? this->num_rows : this->num_columns;
	size_t num_snps = this->transpose ? this->num_columns : this->num_rows;
	sprintf(dest, "%s %7ld %7ld %x %x", matrix_type, num_samples, num_snps, this->individual_file_hash, this->snp_file_hash);
}

int genotype_matrix_write_packed_file(const genotype_matrix *this, const char *filename, bool transpose){
	FILE *output = fopen(filename, "w");

	// first write out header
	int header_size = GENO_HEADER_SIZE;
	if(!transpose){ // header uses row size in packed ancestry map format
		size_t output_row_size = padded_bytes_needed(genotype_matrix_num_samples(this));
		if (header_size < output_row_size){
			header_size = output_row_size;
		}
	}
	//printf("header_size %d\n", header_size);
	char *buffer = malloc(header_size);
	memset(buffer, 0, header_size);
	genotype_matrix_header(this, transpose, buffer);
	size_t bytes_written = fwrite(buffer, sizeof(char), header_size, output);
	if(bytes_written != header_size){
		genotype_matrix_fail("header write failure");
	}
	free(buffer);
	buffer = NULL;

	// if matrix is in same orientation, we simply write it out
	if(this->transpose == transpose){
		bytes_written = fwrite(this->genotypes, sizeof(unsigned char), genotype_matrix_byte_size(this), output);
		if (bytes_written != genotype_matrix_byte_size(this)){
			genotype_matrix_fail("write failure");
		}
	}
	else{ // need to switch rows and columns for output
		size_t output_row_size = padded_bytes_needed(this->num_rows);
		if(!transpose && output_row_size < GENO_HEADER_SIZE){
			output_row_size = GENO_HEADER_SIZE;
		}
		//printf("output row size %d\n", output_row_size);
		buffer = malloc(output_row_size);
		for (size_t column=0; column < this->num_columns; column++){
			//memset(buffer, 0xFF, output_row_size);
			memset(buffer, 0x00, output_row_size);
			for(size_t row = 0; row < this->num_rows; row++){
				bit_array_position read_position;
				row_position(column, &read_position);
				read_position.byte_index += row * this->row_size_bytes;
				unsigned char genotype = genotype_matrix_get_genotype_internal(this, &read_position);
				array_set_genotype((unsigned char *) buffer, row, genotype);
			}
			bytes_written = fwrite(buffer, sizeof(unsigned char), output_row_size, output);
			if(bytes_written != output_row_size){
				genotype_matrix_fail("write failure");
			}
		}
		free(buffer);
	}
	fflush(output);
	fclose(output);
	return 0;
}

void row_position(size_t offset, bit_array_position *position){
	position->byte_index = offset / 4;
	position->bit_shift = 6 - (offset % 4) * 2;
}

// one row is one SNP, all samples
void normal_position(const genotype_matrix *this, size_t sample, size_t snp, bit_array_position *position){
	row_position(sample, position);
	position->byte_index += (this->row_size_bytes * snp);
}

// one row is one sample, all SNPs
void transpose_position(const genotype_matrix *this, size_t sample, size_t snp, bit_array_position *position) {
	row_position(snp, position);
	position->byte_index += (this->row_size_bytes * sample);
}

unsigned char genotype_matrix_get_genotype_internal(const genotype_matrix *this, const bit_array_position *position){
	unsigned char genotype_byte = this->genotypes[position->byte_index];
	return (genotype_byte >> position->bit_shift) & 0x3;
}


int genotype_matrix_get_snp_hash(const genotype_matrix *this){
	return this->snp_file_hash;
}

void genotype_matrix_set_snp_hash(genotype_matrix *this, int hash){
	this->snp_file_hash = hash;
}

int genotype_matrix_get_individual_hash(const genotype_matrix *this){
	return this->individual_file_hash;
}

void genotype_matrix_set_individual_hash(genotype_matrix *this, int hash){
	this->individual_file_hash = hash;
}

size_t genotype_matrix_num_snps(const genotype_matrix *this){
	return (this->transpose) ? this->num_columns : this->num_rows;
}

size_t genotype_matrix_num_samples(const genotype_matrix *this){
	return (this->transpose) ? this->num_rows : this->num_columns;
}

unsigned char genotype_matrix_get_genotype(const genotype_matrix *this, size_t sample, size_t snp){
	bit_array_position position;
	this->get_position(this, sample, snp, &position);

	return genotype_matrix_get_genotype_internal(this, &position);
}

void array_set_genotype_internal(unsigned char *genotypes, const bit_array_position *position, char value){
	if(value < 0){
		value = 0x3;
	}

	unsigned char genotype_byte = genotypes[position->byte_index];
	unsigned char value_mask = 0x3;
	unsigned char existing_mask = 0xFF ^ (value_mask << position->bit_shift);
	genotypes[position->byte_index] = (genotype_byte & existing_mask) | ( (value & value_mask) << position->bit_shift);
}

void array_set_genotype(unsigned char *genotypes, size_t offset, char value){
	bit_array_position position;

	row_position(offset, &position);
	array_set_genotype_internal(genotypes, &position, value);
}

void genotype_matrix_set_genotype(genotype_matrix *this, size_t sample, size_t snp, char value){
	bit_array_position position;
	this->get_position(this, sample, snp, &position);
	array_set_genotype_internal(this->genotypes, &position, value);
}
