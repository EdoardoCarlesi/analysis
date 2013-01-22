#include <stdio.h>

#include "libio/write_io.h"
#include "libcosmo/tinker.h"
#include "libcosmo/mass_function.h"
#include "general_functions.h"
#include "general_variables.h"

int main(int arg, char **argv)
{
	fprintf(stdout, "\nComputing the mass function, starting the main().\n");

	initialize_internal_variables(argv);
	MassFunc.bins  = Settings.n_bins-1; 
	ThMassFunc.bins = Settings.n_bins-1;
	
		// First guess parameters and fits for the different mass functions
	T_mf.A = 0.186; 
	T_mf.a = 1.47; 
	T_mf.b = 2.57; 
	T_mf.c = 1.19;

	Settings.use_one_pk = 1;

		read_halo_file();
		compute_numerical_mass_function();
		print_numerical_mass_function();	

	fprintf(stdout, "\nMass function computed.\n");
	return 0;
}
