#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general_variables.h"
#include "general_functions.h"

#include "libmath/statistics.h"
#include "libio/halo_io.h"

int main(int argc, char **argv){
fprintf(stdout, "\nComputing the Chi^2 distribution of NFW profiles fits for the given halo sample.\n");

	initialize_internal_variables(argv);

	Cosmo.err=0.3;
	Cosmo.OmegaM=1.;

	read_halo_file();

	read_profiles_file();

	calculate_and_sort_chi2();

return 0;
}
