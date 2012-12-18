#include <stdio.h>

#include "libio/halo_io.h"
#include "libio/write_io.h"
#include "libhalo/subhalo_properties.h"
#include "general_variables.h"
#include "general_functions.h"

int main(int argc, char **argv)
{
	fprintf(stderr, "\nmain(). Computing subhalo statistical properties.\n");

		initialize_internal_variables(argv);

		read_halo_file();
	
		compute_subhalo_properties();

		Settings.min_subhaloes=400;

		print_all_subhalo_properties_to_one_file();

	fprintf(stderr, "\nmain(). Subhalo statistical properties computed.\n");

	return 0;
}
