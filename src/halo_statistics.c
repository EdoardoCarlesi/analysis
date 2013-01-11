#include <stdio.h>
#include <stdlib.h>

#include "libio/write_io.h"
#include "libio/halo_io.h"
#include "libhalo/halo_properties.h"

#include "general_functions.h"
#include "general_variables.h"


int main(int argc, char **argv)
{
	fprintf(stderr, "\nmain(). Computing halo statistical properties.\n");

		initialize_internal_variables(argv);

		initialize_halo_properties_structure();

// TODO Check the number of haloes to use!!!!!!!!
		Settings.n_threshold=39;

		read_halo_file();

		compute_halo_properties();

		print_all_halo_properties_to_one_file();
	
	fprintf(stderr, "\nmain(). Halo statistical properties computed.\n");

return 0;
}

