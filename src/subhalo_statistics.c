#include <stdio.h>

#include "libio/io.h"
#include "libhalo/halo.h"

#include "general_def.h"


int main(int argc, char **argv)
{
	fprintf(stdout, "\nmain(). Computing subhalo statistical properties.\n");

		initialize_internal_variables(argv);

		read_halo_file();
	
		compute_subhalo_properties();

		Settings.n_sub_min=400;

		print_all_subhalo_properties_to_one_file();

	fprintf(stdout, "\nmain(). Subhalo statistical properties computed.\n");

	return 0;
}
