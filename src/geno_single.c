#include <stdio.h>
#include <stdlib.h>

#include "geno.h"

// extract a single individual's genotypes from a packed genotype file
int main(int argc, char **argv){
	genotype_matrix geno;
	genotype_matrix selected;

	//printf("reading %s\n", argv[1]);
	size_t ind = strtoull(argv[2], NULL, 10);
	char *output_filename = argv[3];

	genotype_matrix_initialize(&geno);
	int result = genotype_matrix_read_file(&geno, argv[1], 0, NULL, NULL);
	if (result < 0){
		printf("Error reading %s: %d\n", argv[1], result);
		return result;
	}

	size_t num_snps = genotype_matrix_num_snps(&geno);
	size_t num_ind = genotype_matrix_num_samples(&geno);
	if (ind < 0 || ind >= num_ind){
		fprintf(stderr, "Bad index %zu. There are %zu samples.", ind, num_ind);
	}

	genotype_matrix_initialize(&selected);
	bool transpose = true;
	genotype_matrix_allocate(&selected, 1, num_snps, transpose);

	for (size_t snp = 0; snp < num_snps; snp++){
		char value = genotype_matrix_get_genotype(&geno, ind, snp);
		genotype_matrix_set_genotype(&selected, 0, snp, value);
	}
	int snp_hash = genotype_matrix_get_snp_hash(&geno);
	genotype_matrix_set_snp_hash(&selected, snp_hash);

	genotype_matrix_write_packed_file(&selected, output_filename, transpose);

	genotype_matrix_destructor(&geno);
	genotype_matrix_destructor(&selected);
	return 0;
}
