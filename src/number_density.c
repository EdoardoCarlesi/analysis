#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "libio/halo_io.h"
#include "libio/read_io.h"
#include "libio/power_io.h"
#include "libcosmo/number_density.h"
#include "general_variables.h"
#include "general_functions.h"


int main(int argc, char **argv)
{
	fprintf(stderr, "\nmain(). Calculating number density.\n");
	
	initialize_internal_variables(argv);
	
		read_pk_snapshots();
		read_redshift_file();
		get_halo_files_urls();

		compute_number_density();

		//integrate_over_z_range();

	fprintf(stderr, "\nmain(). Number density calculated.\n");

	return 0;
}
