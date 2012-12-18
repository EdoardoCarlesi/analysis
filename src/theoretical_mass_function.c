#include <stdio.h>

#include "libio/power_io.h"
#include "libio/write_io.h"
#include "libcosmo/tinker.h"
#include "libcosmo/cosmological_relations.h"
#include "libcosmo/mass_function.h"

#include "general_variables.h"
#include "general_functions.h"

int main(int arg, char **argv)
{
	fprintf(stderr, "\nComputing the analytical mass function for a given P(k).\n");

	initialize_internal_variables(argv);

	AMF.bins = Settings.n_bins-1;
	
			// First guess parameters and fits for the different mass functions
		T_mf.A = 0.186; T_mf.a = 1.47; T_mf.b = 2.57; T_mf.c = 1.19;

			Settings.rho_0 = default_rho0();

			Settings.use_one_pk = 1;
			init_pks();

			//compute_correlation_function(0.);
			//print_correlation_function();

		compute_theoretical_mass_function(0);

		print_theoretical_mass_function(0.0);

	return 0;
}
