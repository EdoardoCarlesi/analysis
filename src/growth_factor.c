#include <stdio.h>
#include <stdlib.h>

#include "libio/io.h"
#include "libcosmo/cosmo.h"

#include "general_def.h"

int main(int argc, char **argv) 
{
	INFO_MSG("Computing growth factor");

		initialize_internal_variables(argv);

		read_redshift_file();

		init_pks();

		read_pk_snapshots();

		compute_growth_factor();

		print_growth_factor();

	INFO_MSG("Done growth factor computation");

return 0;
}
