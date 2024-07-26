#include <argp.h>
#include <stdio.h>

#include "geno.h"

static const char usage_documentation[] = "Transpose a genotype between a packed ancestry map format and a transpose packed format or vice versa.";

static struct argp_option options[] = {
	{"low-mem",	'm', 0, 0, "Do not allocate memory for the full genotype file. This will substantially reduce memory requirements but will be much slower." },
	{ 0 }
};

static char args_doc[] = "input_genotype_file output_genotype_file";

struct arguments
{
	char *input_file;
	char *output_file;
	int low_mem;
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key){
		case 'm':
			arguments->low_mem = 1;
			break;

		case ARGP_KEY_ARG:
			if (state->arg_num >= 2) // too many arguments
				argp_usage (state);
			else if (state->arg_num == 0)
				arguments->input_file = arg;
			else if (state->arg_num == 1)
				arguments->output_file = arg;
			break;

		case ARGP_KEY_END:
			if(state->arg_num < 2)
				argp_usage(state);
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, usage_documentation };

int main(int argc, char **argv){
	struct arguments arguments;

	// default
	arguments.low_mem = 0;

	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	genotype_matrix geno;

	printf("reading %s\n", arguments.input_file);
	genotype_matrix_initialize(&geno);
	int result;
	if (arguments.low_mem){
		result = genotype_matrix_read_file(&geno, arguments.input_file, 0, NULL, NULL);
	} else{
		result = genotype_matrix_read_file_full(&geno, arguments.input_file, 0, NULL, NULL);
	}
	if (result < 0){
		printf("Error %d", result);
		return result;
	}

	genotype_matrix_write_packed_file(&geno, arguments.output_file, !geno.transpose);

	genotype_matrix_destructor(&geno);
	return 0;
}
