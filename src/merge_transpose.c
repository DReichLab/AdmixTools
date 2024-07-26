#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geno.h"

#define WVERSION "100" 

static const char usage_documentation[] = "Merge tranpose bit-packed genotype files. SNPs for one sample are stored continuously in one row, instead of samples for one SNP being stored in one row as in the packed ancestry format. This is fast because the bytes are already stored in order for each component.";

static struct argp_option options[] = {
	{"input",		'l', "FILE", 0, "File containing list of genotype files or stems (path without .geno) to be merged" },
	{"output",		'o', "FILE", 0, "Output path for genotype file. Include .geno if desired." },
	{"ind_hash",	'i', "int", 0, "Nick's hex hash of individual file" },
	{"snp_hash",	's', "int", 0, "Nick's hex hash of SNP file" },
	{ 0 }
};

struct arguments
{
	char *input_file_list;
	char *output_file;
	int individual_hash;
	int snp_hash;
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key){
		case 'l':
			arguments->input_file_list = arg;
			break;
		case 'o':
			arguments->output_file = arg;
			break;
		case 'i':
			arguments->individual_hash = strtoul(arg, NULL, 16);
			break;
		case 's':
			arguments->snp_hash = strtoul(arg, NULL, 16);
			break;

		case ARGP_KEY_ARG:
			break;

		case ARGP_KEY_END:
			if(arguments->input_file_list == NULL || arguments->output_file == NULL){
				argp_failure(state, 1, 0, "-l and -o are required. See --help");
			}
			//argp_usage (state);
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, 0, usage_documentation };

bool file_exists(const char filename[]){
	FILE *fp = fopen(filename, "r");
	bool exists = (fp != NULL);
	if(fp)
		fclose(fp);
	return exists;
}

// return number of samples added to file
int add_geno(const char filename[], FILE *merged, genotype_matrix *combined_geno){
	genotype_matrix to_add;
	int num_samples = 0;
	// copy filename with enough space for .geno
	int expanded_filename_size = strlen(filename) + 5;
	char *file_to_open = malloc(expanded_filename_size);
	strncpy(file_to_open, filename, expanded_filename_size);
	if (file_to_open[strlen(filename)-1] == '\n')
		file_to_open[strlen(filename)-1] = '\0';

	if (!file_exists(file_to_open)){
		strcat(file_to_open, ".geno");
	}

	printf("reading %s\n", file_to_open);
	genotype_matrix_initialize(&to_add);
	int result = genotype_matrix_read_file(&to_add, file_to_open, 0, NULL, NULL);
	if(result < 0 || !to_add.transpose){
		fprintf(stderr, "Failed to add %s", file_to_open);
		exit(-1);
	}
	if(combined_geno->num_columns <= 0){
		combined_geno->num_columns = to_add.num_columns;
	} else if (combined_geno->num_columns != to_add.num_columns){
		fprintf(stderr, "SNP mismatch %zu %zu %s\n", combined_geno->num_columns, to_add.num_columns, file_to_open);
		exit(-1);
	}
	num_samples = to_add.num_rows;
	fwrite(to_add.genotypes, sizeof(unsigned char), genotype_matrix_byte_size(&to_add), merged);
	genotype_matrix_destructor(&to_add);
	free(file_to_open);
	return num_samples;
}

int main(int argc, char **argv){
	struct arguments arguments;
	arguments.input_file_list = NULL;
	arguments.output_file = NULL;
	arguments.individual_hash = 0;
	arguments.snp_hash = 0;
	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	genotype_matrix combined;
	genotype_matrix_initialize(&combined);
	combined.transpose = true;
	if (arguments.individual_hash)
		combined.individual_file_hash = arguments.individual_hash;
	if (arguments.snp_hash)
		combined.snp_file_hash = arguments.snp_hash;

	FILE *merged;
	merged = fopen(arguments.output_file, "w+");
	if(!merged){
		fprintf(stderr, "failed to open %s\n", arguments.output_file);
		return -1;
	}
	// write header
	char header_buffer[GENO_HEADER_SIZE+1];
	memset(header_buffer, 0, GENO_HEADER_SIZE+1);
	genotype_matrix_header(&combined, true, header_buffer);
	fwrite(header_buffer, sizeof(char), GENO_HEADER_SIZE, merged);

	FILE *list_file;
	size_t len = 2048;
	char *filename = malloc(len);
	ssize_t read;
	list_file = fopen(arguments.input_file_list, "r");
	if(!merged){
		fprintf(stderr, "failed to open %s\n", arguments.input_file_list);
		return -1;
	}
	// combine each component genotype file
	while((read = getline(&filename, &len, list_file)) != -1){
		int samples_added = add_geno(filename, merged, &combined);
		combined.num_rows += samples_added;
	}
	free(filename);

	// go back and rewrite header
	fseek(merged, 0, SEEK_SET);
	genotype_matrix_header(&combined, true, header_buffer);
	fwrite(header_buffer, sizeof(char), GENO_HEADER_SIZE, merged);

	fclose(list_file);
	fclose(merged);

	return 0;
}
